module mod_full3d
  use mod_types
  use mod_config, only: sim_config_type
  use mod_geometry, only: panel_mesh_type
  use mod_operator, only: flatten_current, unflatten_current, vector_norm, apply_full_operator_matrix_free
  implicit none
  private
  public :: run_full_ipo, run_full_ipo_matrix_free

contains

  subroutine run_full_ipo(K,J0,max_order,stop_ratio,orders,ratios,step_times,n_done)
    complex(dp), intent(in) :: K(:,:),J0(:,:,:)
    integer, intent(in) :: max_order
    real(dp), intent(in) :: stop_ratio
    complex(dp), allocatable, intent(out) :: orders(:,:,:,:)
    real(dp), allocatable, intent(out) :: ratios(:),step_times(:)
    integer, intent(out) :: n_done
    complex(dp), allocatable :: previous(:),following(:),cumulative(:)
    real(dp) :: t0,t1
    integer :: r
    allocate(orders(size(J0,1),size(J0,2),size(J0,3),max_order+1))
    allocate(ratios(max_order+1),step_times(max_order+1))
    allocate(previous(size(K,1)),following(size(K,1)),cumulative(size(K,1)))
    orders=(0.0_dp,0.0_dp); ratios=0.0_dp; step_times=0.0_dp
    orders(:,:,:,1)=J0; call flatten_current(J0,previous)
    cumulative=previous; ratios(1)=1.0_dp; n_done=1
    do r=1,max_order
      call cpu_time(t0); following=matmul(K,previous); call cpu_time(t1)
      step_times(r+1)=t1-t0
      cumulative=cumulative+following
      ratios(r+1)=vector_norm(following)/(vector_norm(cumulative)+1.0e-300_dp)
      call unflatten_current(following,orders(:,:,:,r+1)); n_done=r+1
      previous=following
      if (ratios(r+1)<stop_ratio) exit
    end do
    deallocate(previous,following,cumulative)
  end subroutine run_full_ipo

  subroutine run_full_ipo_matrix_free(mesh,cfg,J0,max_order,stop_ratio,orders,ratios,step_times,n_done)
    type(panel_mesh_type), intent(in) :: mesh
    type(sim_config_type), intent(in) :: cfg
    complex(dp), intent(in) :: J0(3,mesh%Q,mesh%M)
    integer, intent(in) :: max_order
    real(dp), intent(in) :: stop_ratio
    complex(dp), allocatable, intent(out) :: orders(:,:,:,:)
    real(dp), allocatable, intent(out) :: ratios(:),step_times(:)
    integer, intent(out) :: n_done
    complex(dp), allocatable :: previous(:,:,:),following(:,:,:),cumulative(:,:,:)
    real(dp) :: t0,t1,numerator,denominator
    integer :: r
    allocate(orders(3,mesh%Q,mesh%M,max_order+1),ratios(max_order+1),step_times(max_order+1))
    allocate(previous(3,mesh%Q,mesh%M),following(3,mesh%Q,mesh%M),cumulative(3,mesh%Q,mesh%M))
    orders=(0.0_dp,0.0_dp); ratios=0.0_dp; step_times=0.0_dp
    previous=J0; cumulative=J0; orders(:,:,:,1)=J0; ratios(1)=1.0_dp; n_done=1
    do r=1,max_order
      call cpu_time(t0); call apply_full_operator_matrix_free(mesh,cfg,previous,following); call cpu_time(t1)
      step_times(r+1)=t1-t0; cumulative=cumulative+following
      numerator=sqrt(sum(abs(following)**2)); denominator=sqrt(sum(abs(cumulative)**2))
      ratios(r+1)=numerator/(denominator+1.0e-300_dp)
      orders(:,:,:,r+1)=following; n_done=r+1; previous=following
      if (ratios(r+1)<stop_ratio) exit
    end do
    deallocate(previous,following,cumulative)
  end subroutine run_full_ipo_matrix_free

end module mod_full3d
