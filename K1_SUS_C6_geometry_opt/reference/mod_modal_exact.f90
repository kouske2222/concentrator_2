module mod_modal_exact
  use mod_types
  use mod_config, only: sim_config_type
  use mod_geometry, only: panel_mesh_type
  use mod_operator_exact, only: pair_map_local
  implicit none
  private
  public :: apply_c6_modal_operator_matrix_free, decompose_current_all_modes, reconstruct_active_modes

contains

  subroutine decompose_current_all_modes(J, Jhat, mode_norms)
    complex(dp), intent(in) :: J(:,:,:)
    complex(dp), intent(out) :: Jhat(3,size(J,2),size(J,3))
    real(dp), intent(out) :: mode_norms(size(J,3))
    complex(dp) :: phase
    integer :: nmode, m, p, q

    nmode = size(J,3)
    Jhat = (0.0_dp, 0.0_dp)
    do m = 0, nmode - 1
      do p = 0, nmode - 1
        phase = exp(-I_C * 2.0_dp * PI * real(m*p,dp) / real(nmode,dp))
        do q = 1, size(J,2)
          Jhat(:,q,m+1) = Jhat(:,q,m+1) + J(:,q,p+1) * phase
        end do
      end do
      mode_norms(m+1) = sqrt(sum(abs(Jhat(:,:,m+1))**2))
    end do
  end subroutine decompose_current_all_modes

  subroutine select_active_modes(mode_norms, mode_policy, mode_tolerance, active_modes)
    real(dp), intent(in) :: mode_norms(:)
    integer, intent(in) :: mode_policy
    real(dp), intent(in) :: mode_tolerance
    integer, allocatable, intent(out) :: active_modes(:)
    logical, allocatable :: active(:)
    integer :: nmode, m, nactive, largest_index
    real(dp) :: largest_mode

    nmode = size(mode_norms)
    allocate(active(nmode))
    active = .false.

    if (mode_policy == 0) then
      active = .true.
    else if (mode_policy == 1) then
      if (nmode >= 2) active(2) = .true.
      if (nmode >= 2) active(nmode) = .true.
    else
      largest_mode = maxval(mode_norms)
      do m = 1, nmode
        active(m) = mode_norms(m) > mode_tolerance * max(largest_mode, 1.0e-300_dp)
      end do
      if (.not. any(active)) then
        largest_index = maxloc(mode_norms, dim=1)
        active(largest_index) = .true.
      end if
    end if

    nactive = count(active)
    allocate(active_modes(nactive))
    nactive = 0
    do m = 0, nmode - 1
      if (active(m+1)) then
        nactive = nactive + 1
        active_modes(nactive) = m
      end if
    end do
    deallocate(active)
  end subroutine select_active_modes

  subroutine gather_active_modes(Jhat, active_modes, Jmode)
    complex(dp), intent(in) :: Jhat(:,:,:)
    integer, intent(in) :: active_modes(:)
    complex(dp), intent(out) :: Jmode(3,size(Jhat,2),size(active_modes))
    integer :: a

    do a = 1, size(active_modes)
      Jmode(:,:,a) = Jhat(:,:,active_modes(a)+1)
    end do
  end subroutine gather_active_modes

  subroutine reconstruct_active_modes(active_modes, Jmode, J)
    integer, intent(in) :: active_modes(:)
    complex(dp), intent(in) :: Jmode(:,:,:)
    complex(dp), intent(out) :: J(:,:,:)
    complex(dp) :: phase
    integer :: nmode, p, q, a, m

    nmode = size(J,3)
    J = (0.0_dp, 0.0_dp)
    do p = 0, nmode - 1
      do a = 1, size(active_modes)
        m = active_modes(a)
        phase = exp(I_C * 2.0_dp * PI * real(m*p,dp) / real(nmode,dp)) / real(nmode,dp)
        do q = 1, size(J,2)
          J(:,q,p+1) = J(:,q,p+1) + Jmode(:,q,a) * phase
        end do
      end do
    end do
  end subroutine reconstruct_active_modes

  subroutine apply_c6_modal_operator_matrix_free(mesh, cfg, active_modes, Jmode, Jmode_next)
    type(panel_mesh_type), intent(in) :: mesh
    type(sim_config_type), intent(in) :: cfg
    integer, intent(in) :: active_modes(:)
    complex(dp), intent(in) :: Jmode(3,mesh%Q,size(active_modes))
    complex(dp), intent(out) :: Jmode_next(3,mesh%Q,size(active_modes))
    complex(dp), allocatable :: phase(:,:)
    complex(dp), allocatable :: acc(:,:)
    complex(dp) :: mapA(3,3), v(3)
    integer :: d, q_target, q_source, a, source_sector

    allocate(phase(mesh%M,size(active_modes)))
    do a = 1, size(active_modes)
      do d = 0, mesh%M - 1
        phase(d+1,a) = exp(I_C * 2.0_dp * PI * real(active_modes(a)*d,dp) / real(mesh%M,dp))
      end do
    end do

    Jmode_next = (0.0_dp, 0.0_dp)
!$omp parallel do schedule(static) default(shared) private(q_target,q_source,d,a,source_sector,mapA,v,acc)
    do q_target = 1, mesh%Q
      allocate(acc(3,size(active_modes)))
      acc = (0.0_dp, 0.0_dp)
      do d = 0, mesh%M - 1
        source_sector = d + 1
        do q_source = 1, mesh%Q
          if (d == 0 .and. q_target == q_source) cycle
          call pair_map_local(mesh, cfg, 1, q_target, source_sector, q_source, mapA)
          if (sum(abs(mapA)) <= tiny(1.0_dp)) cycle
          do a = 1, size(active_modes)
            v = matmul(mapA, Jmode(:,q_source,a)) * phase(d+1,a)
            acc(:,a) = acc(:,a) + v
          end do
        end do
      end do
      do a = 1, size(active_modes)
        Jmode_next(:,q_target,a) = acc(:,a)
      end do
      deallocate(acc)
    end do
!$omp end parallel do

    deallocate(phase)
  end subroutine apply_c6_modal_operator_matrix_free

end module mod_modal_exact
