module mod_operator_exact
  use mod_types
  use mod_config, only: sim_config_type
  use mod_geometry, only: panel_mesh_type, inside_hexagon, local_to_global, global_to_local
  use mod_incident, only: incident_eh_point
  use mod_material, only: surface_current, radiate_panel
  implicit none
  private
  public :: pair_map_local, build_full_operator
  public :: apply_full_operator_matrix_free
  public :: field_from_current, relative_complex_error, vector_norm
  public :: flatten_current, unflatten_current, current_local_to_global

contains

  pure logical function visible_pair(xt,yt,zt,xs,ys,zs,cfg) result(ok)
    real(dp), intent(in) :: xt,yt,zt,xs,ys,zs
    type(sim_config_type), intent(in) :: cfg
    real(dp) :: a,x,y,z
    ! PRECONDITION: endpoints are valid wall centroids from build_c6_mesh.
    ! Both pieces are convex. Their side-plane inequalities are affine
    ! along the segment, so only the common junction needs checking.
    ok=.true.
    z=cfg%base_z+cfg%l_cone
    if (min(zs,zt)>=z .or. max(zs,zt)<=z) return
    a=(z-zs)/(zt-zs)
    x=xs+a*(xt-xs); y=ys+a*(yt-ys)
    ok=inside_hexagon(x,y,z,cfg,0.0_dp)
  end function visible_pair

  subroutine pair_map_local(mesh,cfg,pt,qt,ps,qs,A)
    type(panel_mesh_type), intent(in) :: mesh
    type(sim_config_type), intent(in) :: cfg
    integer, intent(in) :: pt,qt,ps,qs
    complex(dp), intent(out) :: A(3,3)
    real(dp) :: rvec(3),u(3),nt(3),ns(3),R,phi_t,phi_s
    complex(dp) :: coefficient,Jl(3),Jg(3),Hg(3),Jnewg(3),Jnewl(3),Eg(3)
    integer :: column
    A=(0.0_dp,0.0_dp)
    rvec=[mesh%x(pt,qt)-mesh%x(ps,qs),mesh%y(pt,qt)-mesh%y(ps,qs),mesh%z(pt,qt)-mesh%z(ps,qs)]
    ! The self panel is represented by the local PO boundary condition and is
    ! therefore excluded from inter-panel propagation.  No distance softening
    ! or arbitrary near-neighbour deletion is used.
    if (pt==ps .and. qt==qs) return
    R=sqrt(sum(rvec*rvec))
    if (R<=tiny(1.0_dp)) return
    u=rvec/R
    nt=[mesh%nx(pt,qt),mesh%ny(pt,qt),mesh%nz(pt,qt)]
    ns=[mesh%nx(ps,qs),mesh%ny(ps,qs),mesh%nz(ps,qs)]
    if (sum(u*ns)<=1.0e-12_dp .or. sum(u*nt)>=-1.0e-12_dp) return
    if (.not.visible_pair(mesh%x(pt,qt),mesh%y(pt,qt),mesh%z(pt,qt), &
                         mesh%x(ps,qs),mesh%y(ps,qs),mesh%z(ps,qs),cfg)) return
    coefficient=mesh%area(ps,qs)*exp(I_C*cfg%k0*R)/(4.0_dp*PI)
    phi_t=2.0_dp*PI*real(pt-1,dp)/real(mesh%M,dp)
    phi_s=2.0_dp*PI*real(ps-1,dp)/real(mesh%M,dp)
    do column=1,3
      Jl=(0.0_dp,0.0_dp); Jl(column)=(1.0_dp,0.0_dp)
      call local_to_global(phi_s,Jl,Jg)
      call radiate_panel(cfg%k0,cfg%eta0,rvec,mesh%area(ps,qs),ns,Jg,Eg,Hg,coefficient)
      call surface_current(nt,u,Hg,Jnewg)
      call global_to_local(phi_t,Jnewg,Jnewl)
      A(:,column)=Jnewl
    end do
  end subroutine pair_map_local

  subroutine build_full_operator(mesh,cfg,K,elapsed)
    type(panel_mesh_type), intent(in) :: mesh
    type(sim_config_type), intent(in) :: cfg
    complex(dp), allocatable, intent(out) :: K(:,:)
    real(dp), intent(out) :: elapsed
    integer :: pt,qt,ps,qs,it,is
    real(dp) :: t0,t1
    complex(dp) :: A(3,3)
    allocate(K(3*mesh%N,3*mesh%N)); K=(0.0_dp,0.0_dp)
    call cpu_time(t0)
!$omp parallel do collapse(2) schedule(static) default(shared) private(pt,qt,ps,qs,it,is,A)
    do pt=1,mesh%M
      do qt=1,mesh%Q
        it=(pt-1)*mesh%Q+qt
        do ps=1,mesh%M
          do qs=1,mesh%Q
            is=(ps-1)*mesh%Q+qs
            call pair_map_local(mesh,cfg,pt,qt,ps,qs,A)
            K(3*it-2:3*it,3*is-2:3*is)=A
          end do
        end do
      end do
    end do
!$omp end parallel do
    call cpu_time(t1); elapsed=t1-t0
  end subroutine build_full_operator

  subroutine apply_full_operator_matrix_free(mesh,cfg,Jin,Jout)
    type(panel_mesh_type), intent(in) :: mesh
    type(sim_config_type), intent(in) :: cfg
    complex(dp), intent(in) :: Jin(3,mesh%Q,mesh%M)
    complex(dp), intent(out) :: Jout(3,mesh%Q,mesh%M)
    integer :: pt,qt,ps,qs
    complex(dp) :: A(3,3),acc(3)
    Jout=(0.0_dp,0.0_dp)
!$omp parallel do collapse(2) schedule(static) default(shared) private(pt,qt,ps,qs,A,acc)
    do pt=1,mesh%M
      do qt=1,mesh%Q
        acc=(0.0_dp,0.0_dp)
        do ps=1,mesh%M
          do qs=1,mesh%Q
            call pair_map_local(mesh,cfg,pt,qt,ps,qs,A)
            acc=acc+matmul(A,Jin(:,qs,ps))
          end do
        end do
        Jout(:,qt,pt)=acc
      end do
    end do
!$omp end parallel do
  end subroutine apply_full_operator_matrix_free

  subroutine flatten_current(J,flat)
    complex(dp), intent(in) :: J(:,:,:)
    complex(dp), intent(out) :: flat(:)
    integer :: p,q,i
    do p=1,size(J,3); do q=1,size(J,2)
      i=(p-1)*size(J,2)+q; flat(3*i-2:3*i)=J(:,q,p)
    end do; end do
  end subroutine flatten_current

  subroutine unflatten_current(flat,J)
    complex(dp), intent(in) :: flat(:)
    complex(dp), intent(out) :: J(:,:,:)
    integer :: p,q,i
    do p=1,size(J,3); do q=1,size(J,2)
      i=(p-1)*size(J,2)+q; J(:,q,p)=flat(3*i-2:3*i)
    end do; end do
  end subroutine unflatten_current

  real(dp) function vector_norm(v) result(value)
    complex(dp), intent(in) :: v(:)
    value=sqrt(sum(abs(v)**2))
  end function vector_norm

  real(dp) function relative_complex_error(a,b) result(value)
    complex(dp), intent(in) :: a(:),b(:)
    value=vector_norm(a-b)/(vector_norm(b)+1.0e-300_dp)
  end function relative_complex_error

  subroutine current_local_to_global(mesh,Jlocal,Jglobal)
    type(panel_mesh_type), intent(in) :: mesh
    complex(dp), intent(in) :: Jlocal(3,mesh%Q,mesh%M)
    complex(dp), intent(out) :: Jglobal(3,mesh%Q,mesh%M)
    integer :: p,q
    real(dp) :: phi
    do p=0,mesh%M-1
      phi=2.0_dp*PI*real(p,dp)/real(mesh%M,dp)
      do q=1,mesh%Q
        call local_to_global(phi,Jlocal(:,q,p+1),Jglobal(:,q,p+1))
      end do
    end do
  end subroutine current_local_to_global


  subroutine field_from_current(mesh,cfg,Jlocal,xo,yo,zo,E,H,include_incident)
    type(panel_mesh_type), intent(in) :: mesh
    type(sim_config_type), intent(in) :: cfg
    complex(dp), intent(in) :: Jlocal(3,mesh%Q,mesh%M)
    real(dp), intent(in) :: xo(:),yo(:),zo(:)
    complex(dp), intent(out) :: E(3,size(xo)),H(3,size(xo))
    logical, intent(in) :: include_incident
    complex(dp), allocatable :: Jg(:,:,:)
    complex(dp) :: Ei(3),Hi(3)
    real(dp) :: khat(3),rvec(3),ns(3)
    integer :: i,p,q
    allocate(Jg(3,mesh%Q,mesh%M))
    call current_local_to_global(mesh,Jlocal,Jg)
    E=(0.0_dp,0.0_dp); H=(0.0_dp,0.0_dp)
!$omp parallel do default(shared) private(i,p,q,Ei,Hi,khat,rvec,ns)
    do i=1,size(xo)
      if (include_incident) then
        call incident_eh_point(xo(i),yo(i),zo(i),cfg,Ei,Hi,khat)
        E(:,i)=Ei; H(:,i)=Hi
      end if
      do p=1,mesh%M
        do q=1,mesh%Q
          rvec=[xo(i)-mesh%x(p,q),yo(i)-mesh%y(p,q),zo(i)-mesh%z(p,q)]
          ns=[mesh%nx(p,q),mesh%ny(p,q),mesh%nz(p,q)]
          call radiate_panel(cfg%k0,cfg%eta0,rvec,mesh%area(p,q),ns,Jg(:,q,p),Ei,Hi)
          E(:,i)=E(:,i)+Ei; H(:,i)=H(:,i)+Hi
        end do
      end do
    end do
!$omp end parallel do
    deallocate(Jg)
  end subroutine field_from_current
end module mod_operator_exact
