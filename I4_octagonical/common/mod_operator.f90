module mod_operator
  use mod_types
  use mod_config, only: sim_config_type
  use mod_geometry, only: panel_mesh_type, inside_octagon, local_to_global, global_to_local
  use mod_incident, only: incident_eh_point
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
    integer :: k
    real(dp) :: a,x,y,z
    ok=.true.
    do k=1,cfg%visibility_samples
      a=real(k,dp)/real(cfg%visibility_samples+1,dp)
      x=xs+a*(xt-xs); y=ys+a*(yt-ys); z=zs+a*(zt-zs)
      if (.not.inside_octagon(x,y,z,cfg,0.0_dp)) then
        ok=.false.; return
      end if
    end do
  end function visible_pair

  subroutine pair_map_local(mesh,cfg,pt,qt,ps,qs,A)
    type(panel_mesh_type), intent(in) :: mesh
    type(sim_config_type), intent(in) :: cfg
    integer, intent(in) :: pt,qt,ps,qs
    complex(dp), intent(out) :: A(3,3)
    real(dp) :: rvec(3),u(3),nt(3),ns(3),R,phi_t,phi_s
    complex(dp) :: coefficient,Jl(3),Jg(3),Hg(3),Jnewg(3),Jnewl(3)
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
    coefficient=mesh%area(ps,qs)*exp(I_C*cfg%k0*R)*(I_C*cfg%k0/R-1.0_dp/(R*R))/(4.0_dp*PI)
    phi_t=2.0_dp*PI*real(pt-1,dp)/real(mesh%M,dp)
    phi_s=2.0_dp*PI*real(ps-1,dp)/real(mesh%M,dp)
    do column=1,3
      Jl=(0.0_dp,0.0_dp); Jl(column)=(1.0_dp,0.0_dp)
      call local_to_global(phi_s,Jl,Jg)
      Hg(1)=coefficient*(u(2)*Jg(3)-u(3)*Jg(2))
      Hg(2)=coefficient*(u(3)*Jg(1)-u(1)*Jg(3))
      Hg(3)=coefficient*(u(1)*Jg(2)-u(2)*Jg(1))
      Jnewg(1)=2.0_dp*(nt(2)*Hg(3)-nt(3)*Hg(2))
      Jnewg(2)=2.0_dp*(nt(3)*Hg(1)-nt(1)*Hg(3))
      Jnewg(3)=2.0_dp*(nt(1)*Hg(2)-nt(2)*Hg(1))
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
    complex(dp) :: Ei(3),Hi(3),phase,a,b,ch,j(3),ju
    real(dp) :: khat(3),rvec(3),u(3),R,invR
    integer :: i,p,q
    allocate(Jg(3,mesh%Q,mesh%M)); call current_local_to_global(mesh,Jlocal,Jg)
    E=(0.0_dp,0.0_dp); H=(0.0_dp,0.0_dp)
    do i=1,size(xo)
      if (include_incident) then
        call incident_eh_point(xo(i),yo(i),zo(i),cfg,Ei,Hi,khat)
        E(:,i)=Ei; H(:,i)=Hi
      end if
      do p=1,mesh%M; do q=1,mesh%Q
        rvec=[xo(i)-mesh%x(p,q),yo(i)-mesh%y(p,q),zo(i)-mesh%z(p,q)]
        R=sqrt(sum(rvec*rvec))
        ! Observation points are masked away from the wall by the caller.
        ! This guard protects only the exact zero-distance case.
        if (R<=tiny(1.0_dp)) cycle
        invR=1.0_dp/R; u=rvec*invR; j=Jg(:,q,p); ju=sum(j*cmplx(u,0.0_dp,kind=dp))
        phase=exp(I_C*cfg%k0*R)*mesh%area(p,q)
        a=(-I_C*cfg%k0)*invR
        b=invR*invR+I_C*invR*invR*invR/cfg%k0
        E(:,i)=E(:,i)+cfg%eta0/(4.0_dp*PI)*phase*(a*(u*ju-j)+b*(3.0_dp*u*ju-j))
        ch=(I_C*cfg%k0*invR-invR*invR)
        H(1,i)=H(1,i)+phase*ch*(u(2)*j(3)-u(3)*j(2))/(4.0_dp*PI)
        H(2,i)=H(2,i)+phase*ch*(u(3)*j(1)-u(1)*j(3))/(4.0_dp*PI)
        H(3,i)=H(3,i)+phase*ch*(u(1)*j(2)-u(2)*j(1))/(4.0_dp*PI)
      end do; end do
    end do
    deallocate(Jg)
  end subroutine field_from_current

end module mod_operator
