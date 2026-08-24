module mod_c6_modal
  use mod_types
  use mod_config, only: sim_config_type
  use mod_geometry, only: panel_mesh_type
  use mod_operator, only: pair_map_local
  implicit none
  private
  public :: run_c6_modal_ipo, decompose_current_all_modes, reconstruct_active_modes

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

  subroutine checkpoint_meta_filename(output_dir, filename)
    character(len=*), intent(in) :: output_dir
    character(len=512), intent(out) :: filename
    write(filename,'(A,"/c6_matrix_free_checkpoint_meta.bin")') trim(output_dir)
  end subroutine checkpoint_meta_filename

  subroutine checkpoint_order_filename(output_dir, order, filename)
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: order
    character(len=512), intent(out) :: filename
    write(filename,'(A,"/c6_matrix_free_order_",I6.6,".bin")') trim(output_dir), order
  end subroutine checkpoint_order_filename

  subroutine save_checkpoint_meta(output_dir, q_count, sector_count, max_order, last_order, &
                                  active_modes, ratios, step_times)
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: q_count, sector_count, max_order, last_order
    integer, intent(in) :: active_modes(:)
    real(dp), intent(in) :: ratios(:), step_times(:)
    character(len=8) :: magic
    character(len=512) :: filename
    integer :: u

    magic = 'C6MF002 '
    call checkpoint_meta_filename(output_dir, filename)
    open(newunit=u, file=trim(filename), status='replace', form='unformatted', action='write')
    write(u) magic
    write(u) q_count
    write(u) sector_count
    write(u) max_order
    write(u) last_order
    write(u) size(active_modes)
    write(u) active_modes
    write(u) ratios(1:max_order+1)
    write(u) step_times(1:max_order+1)
    close(u)
  end subroutine save_checkpoint_meta

  subroutine load_checkpoint_meta(output_dir, q_expected, sector_expected, max_order_expected, &
                                  active_modes_expected, last_order, ratios, step_times, ok)
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: q_expected, sector_expected, max_order_expected
    integer, intent(in) :: active_modes_expected(:)
    integer, intent(out) :: last_order
    real(dp), intent(out) :: ratios(:), step_times(:)
    logical, intent(out) :: ok
    character(len=8) :: magic
    character(len=512) :: filename
    integer, allocatable :: active_modes_file(:)
    integer :: u, q_file, sector_file, max_order_file, nactive_file

    call checkpoint_meta_filename(output_dir, filename)
    inquire(file=trim(filename), exist=ok)
    if (.not. ok) then
      last_order = -1
      return
    end if

    open(newunit=u, file=trim(filename), status='old', form='unformatted', action='read')
    read(u) magic
    read(u) q_file
    read(u) sector_file
    read(u) max_order_file
    read(u) last_order
    read(u) nactive_file
    allocate(active_modes_file(nactive_file))
    read(u) active_modes_file
    if (magic /= 'C6MF002 ') stop 'Invalid C6 matrix-free checkpoint meta magic.'
    if (q_file /= q_expected) stop 'C6 matrix-free checkpoint panel count mismatch.'
    if (sector_file /= sector_expected) stop 'C6 matrix-free checkpoint sector count mismatch.'
    if (max_order_file /= max_order_expected) stop 'C6 matrix-free checkpoint max_order mismatch.'
    if (last_order < 0 .or. last_order > max_order_expected) stop 'Invalid C6 matrix-free checkpoint last_order.'
    if (nactive_file /= size(active_modes_expected)) stop 'C6 matrix-free checkpoint active-mode count mismatch.'
    if (any(active_modes_file /= active_modes_expected)) stop 'C6 matrix-free checkpoint active-mode list mismatch.'
    read(u) ratios(1:max_order_expected+1)
    read(u) step_times(1:max_order_expected+1)
    close(u)
    deallocate(active_modes_file)
  end subroutine load_checkpoint_meta

  subroutine save_checkpoint(output_dir, order, active_modes, Jmode)
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: order
    integer, intent(in) :: active_modes(:)
    complex(dp), intent(in) :: Jmode(:,:,:)
    character(len=8) :: magic
    character(len=512) :: filename
    integer :: u

    magic = 'C6MFO02 '
    call checkpoint_order_filename(output_dir, order, filename)
    open(newunit=u, file=trim(filename), status='replace', form='unformatted', action='write')
    write(u) magic
    write(u) order
    write(u) size(Jmode,2)
    write(u) size(active_modes)
    write(u) active_modes
    write(u) Jmode
    close(u)
  end subroutine save_checkpoint

  subroutine load_checkpoint(output_dir, order, active_modes, Jmode)
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: order
    integer, intent(in) :: active_modes(:)
    complex(dp), intent(out) :: Jmode(:,:,:)
    character(len=8) :: magic
    character(len=512) :: filename
    integer, allocatable :: active_modes_file(:)
    integer :: u, order_file, q_file, nactive_file

    call checkpoint_order_filename(output_dir, order, filename)
    open(newunit=u, file=trim(filename), status='old', form='unformatted', action='read')
    read(u) magic
    read(u) order_file
    read(u) q_file
    read(u) nactive_file
    allocate(active_modes_file(nactive_file))
    read(u) active_modes_file
    if (magic /= 'C6MFO02 ') stop 'Invalid C6 matrix-free checkpoint order magic.'
    if (order_file /= order) stop 'C6 matrix-free checkpoint order mismatch.'
    if (q_file /= size(Jmode,2)) stop 'C6 matrix-free checkpoint order panel count mismatch.'
    if (nactive_file /= size(active_modes)) stop 'C6 matrix-free checkpoint order active-mode count mismatch.'
    if (any(active_modes_file /= active_modes)) stop 'C6 matrix-free checkpoint order active-mode list mismatch.'
    read(u) Jmode
    close(u)
    deallocate(active_modes_file)
  end subroutine load_checkpoint

  subroutine run_c6_modal_ipo(mesh, cfg, J0, max_order, stop_ratio, orders, ratios, step_times, &
                              initial_mode_norms, final_mode_norms, n_done)
    type(panel_mesh_type), intent(in) :: mesh
    type(sim_config_type), intent(in) :: cfg
    complex(dp), intent(in) :: J0(3,mesh%Q,mesh%M)
    integer, intent(in) :: max_order
    real(dp), intent(in) :: stop_ratio
    complex(dp), allocatable, intent(out) :: orders(:,:,:,:)
    real(dp), allocatable, intent(out) :: ratios(:), step_times(:)
    real(dp), intent(out) :: initial_mode_norms(mesh%M), final_mode_norms(mesh%M)
    integer, intent(out) :: n_done
    complex(dp), allocatable :: Jhat_all(:,:,:), Jmode(:,:,:), Jmode_next(:,:,:)
    complex(dp), allocatable :: Jmode_saved(:,:,:), cumulative_mode(:,:,:)
    integer, allocatable :: active_modes(:)
    real(dp) :: t0, t1, numerator, denominator
    integer :: order, loaded_last, a
    logical :: have_checkpoint

    allocate(Jhat_all(3,mesh%Q,mesh%M))
    call decompose_current_all_modes(J0, Jhat_all, initial_mode_norms)
    call select_active_modes(initial_mode_norms, cfg%mode_policy, cfg%mode_tolerance, active_modes)

    allocate(orders(3,mesh%Q,mesh%M,max_order+1))
    allocate(ratios(max_order+1), step_times(max_order+1))
    allocate(Jmode(3,mesh%Q,size(active_modes)), Jmode_next(3,mesh%Q,size(active_modes)), &
             Jmode_saved(3,mesh%Q,size(active_modes)), cumulative_mode(3,mesh%Q,size(active_modes)))

    call gather_active_modes(Jhat_all, active_modes, Jmode)
    orders = (0.0_dp, 0.0_dp)
    ratios = 0.0_dp
    step_times = 0.0_dp
    final_mode_norms = 0.0_dp

    call load_checkpoint_meta(trim(cfg%output_dir), mesh%Q, mesh%M, max_order, active_modes, loaded_last, &
                              ratios, step_times, have_checkpoint)

    if (have_checkpoint) then
      cumulative_mode = (0.0_dp, 0.0_dp)
      do order = 0, loaded_last
        call load_checkpoint(trim(cfg%output_dir), order, active_modes, Jmode_saved)
        cumulative_mode = cumulative_mode + Jmode_saved
        call reconstruct_active_modes(active_modes, Jmode_saved, orders(:,:,:,order+1))
      end do
      Jmode = Jmode_saved
      n_done = loaded_last + 1
    else
      call reconstruct_active_modes(active_modes, Jmode, orders(:,:,:,1))
      cumulative_mode = Jmode
      ratios(1) = 1.0_dp
      step_times(1) = 0.0_dp
      n_done = 1
      call save_checkpoint(trim(cfg%output_dir), 0, active_modes, Jmode)
      call save_checkpoint_meta(trim(cfg%output_dir), mesh%Q, mesh%M, max_order, 0, active_modes, ratios, step_times)
      loaded_last = 0
    end if

    if (loaded_last < max_order) then
      do order = loaded_last + 1, max_order
        call cpu_time(t0)
        call apply_c6_modal_operator_matrix_free(mesh, cfg, active_modes, Jmode, Jmode_next)
        call cpu_time(t1)
        step_times(order+1) = t1 - t0
        call reconstruct_active_modes(active_modes, Jmode_next, orders(:,:,:,order+1))
        cumulative_mode = cumulative_mode + Jmode_next
        numerator = sqrt(sum(abs(Jmode_next)**2))
        denominator = sqrt(sum(abs(cumulative_mode)**2))
        ratios(order+1) = numerator / (denominator + 1.0e-300_dp)
        n_done = order + 1
        Jmode = Jmode_next
        call save_checkpoint(trim(cfg%output_dir), order, active_modes, Jmode)
        call save_checkpoint_meta(trim(cfg%output_dir), mesh%Q, mesh%M, max_order, order, active_modes, ratios, step_times)
        if (ratios(order+1) < stop_ratio) exit
      end do
    end if

    final_mode_norms = 0.0_dp
    do a = 1, size(active_modes)
      final_mode_norms(active_modes(a)+1) = sqrt(sum(abs(Jmode(:,:,a))**2))
    end do

    deallocate(Jhat_all, Jmode, Jmode_next, Jmode_saved, cumulative_mode, active_modes)
  end subroutine run_c6_modal_ipo

end module mod_c6_modal
