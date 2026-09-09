! One immutable geometry/material context per process. Not reentrant.
module mod_bridge
  use iso_c_binding
  use mod_types
  use mod_config
  use mod_geometry
  use mod_material
  use mod_operator, only: pair_map_local
  use mod_c6_modal
  implicit none
  type(sim_config_type), save :: cfg
  type(panel_mesh_type), save :: mesh
contains
  subroutine hm_near(d,nref,nblock,ncol,rowptr,refs,colptr,cols,x,y) bind(C)
    integer(c_int), value :: d,nref,nblock,ncol
    integer(c_int), intent(in) :: rowptr(mesh%Q+1),refs(nref),colptr(nblock+1),cols(ncol)
    complex(c_double_complex), intent(in) :: x(3,mesh%Q,6)
    complex(c_double_complex), intent(out) :: y(3,mesh%Q,6)
    complex(dp) :: a(3,3),acc(3,6)
    integer :: i,j,ref,b,m
    ! One thread owns each target: no atomic sums, no Q^2 adjacency table.
!$omp parallel do default(shared) private(i,j,ref,b,m,a,acc)
    do i=1,mesh%Q
      acc=(0.0_dp,0.0_dp)
      do ref=rowptr(i)+1,rowptr(i+1)
        b=refs(ref)+1
        do j=colptr(b)+1,colptr(b+1)
          call pair_map_local(mesh,cfg,1,i,d+1,cols(j)+1,a)
          if(sum(abs(a))<=tiny(1.0_dp)) cycle
          do m=1,6
            acc(:,m)=acc(:,m)+matmul(a,x(:,cols(j)+1,m))
          end do
        end do
      end do
      y(:,i,:)=acc
    end do
!$omp end parallel do
  end subroutine

  subroutine hm_panel_matrix(d,nr,nc,rows,cols,aout) bind(C)
    integer(c_int), value :: d,nr,nc
    integer(c_int), intent(in) :: rows(nr),cols(nc)
    complex(c_double_complex), intent(out) :: aout(3,nc,3,nr)
    complex(dp) :: a(3,3)
    integer :: i,j,c,t
!$omp parallel do default(shared) private(i,j,c,t,a) if(nr*nc>256)
    do i=1,nr
      do j=1,nc
        call pair_map_local(mesh,cfg,1,rows(i)+1,d+1,cols(j)+1,a)
        do t=1,3
          do c=1,3
            aout(c,j,t,i)=a(t,c)
          end do
        end do
      end do
    end do
!$omp end parallel do
  end subroutine
  subroutine hm_apply_direct(j,jout) bind(C)
    complex(c_double_complex), intent(in) :: j(3,mesh%Q,6)
    complex(c_double_complex), intent(out) :: jout(3,mesh%Q,6)
    complex(dp), allocatable :: modes(:,:,:),next(:,:,:)
    real(dp) :: norms(6)
    allocate(modes(3,mesh%Q,6),next(3,mesh%Q,6))
    call decompose_current_all_modes(j,modes,norms)
    call apply_c6_modal_operator_matrix_free(mesh,cfg,[0,1,2,3,4,5],modes,next)
    call reconstruct_active_modes([0,1,2,3,4,5],next,jout)
  end subroutine
  subroutine hm_init(filename,nchar,q,qcone,k) bind(C)
    integer(c_int), value :: nchar
    character(c_char), intent(in) :: filename(nchar)
    integer(c_int), intent(out) :: q,qcone
    real(c_double), intent(out) :: k
    character(len=:), allocatable :: path
    integer :: i
    allocate(character(len=nchar)::path)
    do i=1,nchar
      path(i:i)=filename(i)
    end do
    call load_config(path,cfg)
    call load_material(path,cfg)
    call build_c6_mesh(cfg,mesh)
    q=mesh%Q; qcone=2*cfg%n_face*cfg%n_z_cone; k=cfg%k0
  end subroutine
  subroutine hm_coords(xyz) bind(C)
    real(c_double), intent(out) :: xyz(3,mesh%Q,6)
    integer :: p,i
    do p=1,6
      do i=1,mesh%Q
        xyz(:,i,p)=[mesh%x(p,i),mesh%y(p,i),mesh%z(p,i)]
      end do
    end do
  end subroutine
  subroutine hm_entries(d,n,rows,cols,values) bind(C)
    integer(c_int), value :: d,n
    integer(c_int), intent(in) :: rows(n),cols(n)
    complex(c_double_complex), intent(out) :: values(n)
    complex(dp) :: a(3,3)
    integer :: i,r,c
!$omp parallel do default(shared) private(i,r,c,a) if(n>512)
    do i=1,n
      r=rows(i); c=cols(i)
      call pair_map_local(mesh,cfg,1,r/3+1,d+1,c/3+1,a)
      values(i)=a(mod(r,3)+1,mod(c,3)+1)
    end do
!$omp end parallel do
  end subroutine
  subroutine hm_direct(d,nr,nc,rows,cols,x,y) bind(C)
    integer(c_int), value :: d,nr,nc
    integer(c_int), intent(in) :: rows(nr),cols(nc)
    complex(c_double_complex), intent(in) :: x(3,nc,6)
    complex(c_double_complex), intent(out) :: y(3,nr,6)
    complex(dp) :: a(3,3),acc(3,6)
    integer :: i,j,m
!$omp parallel do default(shared) private(i,j,m,a,acc) if(nr*nc>256)
    do i=1,nr
      acc=(0.0_dp,0.0_dp)
      do j=1,nc
        call pair_map_local(mesh,cfg,1,rows(i)+1,d+1,cols(j)+1,a)
        do m=1,6
          acc(:,m)=acc(:,m)+matmul(a,x(:,j,m))
        end do
      end do
      y(:,i,:)=acc
    end do
!$omp end parallel do
  end subroutine
end module
