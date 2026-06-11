module mod_currents_solver
  use mod_types
  use mod_checkpoint, only: ensure_checkpoint_dir, save_ipo_checkpoint_meta, &
       load_ipo_checkpoint_meta, save_ipo_checkpoint_order, load_ipo_checkpoint_order
  use mod_incident, only: calc_gaussian_incident_eh_points, calc_surface_currents_from_incident_eh
  use mod_field_integrals, only: calc_eh_from_surface_currents_to_points
  implicit none
  private
  public :: propagate_currents_to_new_surface_order
  public :: build_iterative_po_currents_fast
  public :: sum_currents_by_order

contains

  subroutine propagate_currents_to_new_surface_order( &
      src_Jx, src_Jy, src_Jz, &
      src_x, src_y, src_z, src_dS, src_panel_len, &
      tgt_x, tgt_y, tgt_z, tgt_nx, tgt_ny, tgt_nz, &
      k0, eta0, exclusion_factor, apply_illumination_mask, &
      Jx_new, Jy_new, Jz_new, &
      current_order, total_orders, show_progress, n_progress_updates)
    complex(dp), intent(in) :: src_Jx(:), src_Jy(:), src_Jz(:)
    real(dp), intent(in) :: src_x(:), src_y(:), src_z(:), src_dS(:), src_panel_len(:)
    real(dp), intent(in) :: tgt_x(:), tgt_y(:), tgt_z(:), tgt_nx(:), tgt_ny(:), tgt_nz(:)
    real(dp), intent(in) :: k0, eta0, exclusion_factor
    logical, intent(in) :: apply_illumination_mask
    complex(dp), intent(out) :: Jx_new(:), Jy_new(:), Jz_new(:)
    integer, intent(in), optional :: current_order, total_orders, n_progress_updates
    logical, intent(in), optional :: show_progress

    integer :: N_tgt
    integer :: i, ibeg, iend, block_size, n_updates
    integer :: percent_now, percent_prev
    complex(dp), allocatable :: Ex(:), Ey(:), Ez(:), Hx(:), Hy(:), Hz(:)
    real(dp) :: sx, sy, sz
    logical :: illuminated, do_show

    N_tgt = size(tgt_x)

    do_show = .false.
    if (present(show_progress)) do_show = show_progress

    n_updates = 20
    if (present(n_progress_updates)) n_updates = max(1, n_progress_updates)

    if (do_show) then
      block_size = max(1, ceiling(real(N_tgt, dp) / real(n_updates, dp)))
    else
      block_size = N_tgt
    end if

    allocate(Ex(N_tgt), Ey(N_tgt), Ez(N_tgt), Hx(N_tgt), Hy(N_tgt), Hz(N_tgt))

    if (.not. do_show) then
      call calc_eh_from_surface_currents_to_points( &
          src_Jx, src_Jy, src_Jz, &
          src_x, src_y, src_z, src_dS, src_panel_len, &
          tgt_x, tgt_y, tgt_z, &
          k0, eta0, exclusion_factor, &
          Ex, Ey, Ez, Hx, Hy, Hz )
    else
      percent_prev = -1

      do ibeg = 1, N_tgt, block_size
        iend = min(ibeg + block_size - 1, N_tgt)

        call calc_eh_from_surface_currents_to_points( &
            src_Jx, src_Jy, src_Jz, &
            src_x, src_y, src_z, src_dS, src_panel_len, &
            tgt_x(ibeg:iend), tgt_y(ibeg:iend), tgt_z(ibeg:iend), &
            k0, eta0, exclusion_factor, &
            Ex(ibeg:iend), Ey(ibeg:iend), Ez(ibeg:iend), &
            Hx(ibeg:iend), Hy(ibeg:iend), Hz(ibeg:iend) )

        percent_now = int(100.0_dp * real(iend, dp) / real(N_tgt, dp))
        if (percent_now /= percent_prev) then
          if (present(current_order) .and. present(total_orders)) then
            write(*,'(A,I0,A,I0,A,I3,A)') &
                '    IPO current: order ', current_order, '/', total_orders, '  ', percent_now, '%'
          else
            write(*,'(A,I3,A)') '    IPO current: ', percent_now, '%'
          end if
          percent_prev = percent_now
        end if
      end do
    end if

    Jx_new = (0.0_dp, 0.0_dp)
    Jy_new = (0.0_dp, 0.0_dp)
    Jz_new = (0.0_dp, 0.0_dp)

!$omp parallel do schedule(static) default(shared) private(i,sx,sy,sz,illuminated)
    do i = 1, N_tgt
      illuminated = .true.

      if (apply_illumination_mask) then
        sx = 0.5_dp * real(Ey(i)*conjg(Hz(i)) - Ez(i)*conjg(Hy(i)), dp)
        sy = 0.5_dp * real(Ez(i)*conjg(Hx(i)) - Ex(i)*conjg(Hz(i)), dp)
        sz = 0.5_dp * real(Ex(i)*conjg(Hy(i)) - Ey(i)*conjg(Hx(i)), dp)
        illuminated = (sx*tgt_nx(i) + sy*tgt_ny(i) + sz*tgt_nz(i)) < 0.0_dp
      end if

      if (illuminated) then
        Jx_new(i) = 2.0_dp * (tgt_ny(i)*Hz(i) - tgt_nz(i)*Hy(i))
        Jy_new(i) = 2.0_dp * (tgt_nz(i)*Hx(i) - tgt_nx(i)*Hz(i))
        Jz_new(i) = 2.0_dp * (tgt_nx(i)*Hy(i) - tgt_ny(i)*Hx(i))
      end if
    end do
!$omp end parallel do

    deallocate(Ex, Ey, Ez, Hx, Hy, Hz)
  end subroutine propagate_currents_to_new_surface_order

  subroutine build_iterative_po_currents_fast( &
      wall_x, wall_y, wall_z, wall_nx, wall_ny, wall_nz, wall_dS, wall_panel_len, &
      max_reflection_order, k0, eta0, w0_src, z0_src, zR_src, exclusion_factor_wall, stop_ratio, &
      currents_x, currents_y, currents_z, ratios, order_build_times, n_orders_done, &
      show_progress, n_progress_updates, checkpoint_dir, restart_from_checkpoint)
    real(dp), intent(in) :: wall_x(:), wall_y(:), wall_z(:), wall_nx(:), wall_ny(:), wall_nz(:)
    real(dp), intent(in) :: wall_dS(:), wall_panel_len(:)
    integer, intent(in) :: max_reflection_order
    real(dp), intent(in) :: k0, eta0, w0_src, z0_src, zR_src, exclusion_factor_wall, stop_ratio
    complex(dp), allocatable, intent(out) :: currents_x(:,:), currents_y(:,:), currents_z(:,:)
    real(dp), allocatable, intent(out) :: ratios(:), order_build_times(:)
    integer, intent(out) :: n_orders_done
    logical, intent(in), optional :: show_progress
    integer, intent(in), optional :: n_progress_updates
    character(len=*), intent(in), optional :: checkpoint_dir
    logical, intent(in), optional :: restart_from_checkpoint

    integer :: N, order, start_order, last_order_loaded
    integer :: n_updates
    logical :: do_show, use_checkpoint, use_restart, meta_ok
    character(len=512) :: ckpt_dir
    complex(dp), allocatable :: Ex0(:), Ey0(:), Ez0(:), Hx0(:), Hy0(:), Hz0(:)
    complex(dp), allocatable :: Jx_new(:), Jy_new(:), Jz_new(:)
    real(dp) :: rms1, rms_now, t_order0, t_order1

    do_show = .false.
    if (present(show_progress)) do_show = show_progress

    n_updates = 20
    if (present(n_progress_updates)) n_updates = max(1, n_progress_updates)

    use_checkpoint = .false.
    ckpt_dir = ''
    if (present(checkpoint_dir)) then
      if (len_trim(checkpoint_dir) > 0) then
        use_checkpoint = .true.
        ckpt_dir = trim(checkpoint_dir)
        call ensure_checkpoint_dir(ckpt_dir)
      end if
    end if

    use_restart = .false.
    if (present(restart_from_checkpoint)) use_restart = restart_from_checkpoint

    N = size(wall_x)

    allocate(currents_x(N, max_reflection_order))
    allocate(currents_y(N, max_reflection_order))
    allocate(currents_z(N, max_reflection_order))
    allocate(ratios(max_reflection_order))
    allocate(order_build_times(max_reflection_order))

    currents_x = (0.0_dp, 0.0_dp)
    currents_y = (0.0_dp, 0.0_dp)
    currents_z = (0.0_dp, 0.0_dp)
    ratios = 0.0_dp
    order_build_times = 0.0_dp

    n_orders_done = 0
    start_order = 1
    last_order_loaded = 0
    meta_ok = .false.

    if (use_checkpoint .and. use_restart) then
      call load_ipo_checkpoint_meta( &
          ckpt_dir, N, max_reflection_order, last_order_loaded, ratios, order_build_times, meta_ok)

      if (meta_ok .and. last_order_loaded >= 1) then
        do order = 1, last_order_loaded
          call load_ipo_checkpoint_order( &
              ckpt_dir, order, currents_x(:,order), currents_y(:,order), currents_z(:,order))
        end do

        n_orders_done = last_order_loaded
        start_order = last_order_loaded + 1

        rms1 = sqrt(sum(abs(currents_x(:,1))**2 + abs(currents_y(:,1))**2 + abs(currents_z(:,1))**2) / real(N, dp))

        if (do_show) then
          write(*,'(A,I0)') '  checkpoint restart: last_order_loaded = ', last_order_loaded
        end if
      end if
    end if

    if (n_orders_done >= max_reflection_order) then
      if (do_show) write(*,'(A)') '  all IPO orders already available in checkpoint'
      return
    end if

    if (n_orders_done > 0) then
      if (stop_ratio > 0.0_dp .and. ratios(n_orders_done) < stop_ratio) then
        if (do_show) write(*,'(A)') '  checkpoint already reached stop_ratio'
        return
      end if
    end if

    if (start_order <= 1) then
      allocate(Ex0(N), Ey0(N), Ez0(N), Hx0(N), Hy0(N), Hz0(N))

      if (do_show) write(*,'(A,I0,A,I0,A)') '  order ', 1, '/', max_reflection_order, ' start'

      call cpu_time(t_order0)
      call calc_gaussian_incident_eh_points( &
          wall_x, wall_y, wall_z, w0_src, z0_src, k0, zR_src, eta0, Ex0, Ey0, Ez0, Hx0, Hy0, Hz0)
      call calc_surface_currents_from_incident_eh( &
          wall_nx, wall_ny, wall_nz, Ex0, Ey0, Ez0, Hx0, Hy0, Hz0, &
          .true., currents_x(:,1), currents_y(:,1), currents_z(:,1))
      call cpu_time(t_order1)

      order_build_times(1) = t_order1 - t_order0
      rms1 = sqrt(sum(abs(currents_x(:,1))**2 + abs(currents_y(:,1))**2 + abs(currents_z(:,1))**2) / real(N, dp))
      ratios(1) = 1.0_dp
      n_orders_done = 1

      if (use_checkpoint) then
        call save_ipo_checkpoint_order(ckpt_dir, 1, currents_x(:,1), currents_y(:,1), currents_z(:,1))
        call save_ipo_checkpoint_meta(ckpt_dir, N, max_reflection_order, n_orders_done, ratios, order_build_times)
      end if

      if (do_show) then
        write(*,'(A,I0,A,I0,A,ES12.4)') '  order ', 1, '/', max_reflection_order, ' done, ratio = ', ratios(1)
      end if

      deallocate(Ex0, Ey0, Ez0, Hx0, Hy0, Hz0)
      start_order = 2
    end if

    if (n_orders_done >= max_reflection_order) return
    if (stop_ratio > 0.0_dp .and. ratios(n_orders_done) < stop_ratio) return

    allocate(Jx_new(N), Jy_new(N), Jz_new(N))

    do order = start_order, max_reflection_order
      if (do_show) write(*,'(A,I0,A,I0,A)') '  order ', order, '/', max_reflection_order, ' start'

      call cpu_time(t_order0)
      call propagate_currents_to_new_surface_order( &
          currents_x(:,order-1), currents_y(:,order-1), currents_z(:,order-1), &
          wall_x, wall_y, wall_z, wall_dS, wall_panel_len, &
          wall_x, wall_y, wall_z, wall_nx, wall_ny, wall_nz, &
          k0, eta0, exclusion_factor_wall, .true., &
          Jx_new, Jy_new, Jz_new, &
          current_order=order, total_orders=max_reflection_order, &
          show_progress=do_show, n_progress_updates=n_updates)
      call cpu_time(t_order1)

      order_build_times(order) = t_order1 - t_order0
      currents_x(:,order) = Jx_new
      currents_y(:,order) = Jy_new
      currents_z(:,order) = Jz_new

      rms_now = sqrt(sum(abs(Jx_new)**2 + abs(Jy_new)**2 + abs(Jz_new)**2) / real(N, dp))
      ratios(order) = rms_now / (rms1 + 1.0e-30_dp)
      n_orders_done = order

      if (use_checkpoint) then
        call save_ipo_checkpoint_order(ckpt_dir, order, currents_x(:,order), currents_y(:,order), currents_z(:,order))
        call save_ipo_checkpoint_meta(ckpt_dir, N, max_reflection_order, n_orders_done, ratios, order_build_times)
      end if

      if (do_show) then
        write(*,'(A,I0,A,I0,A,ES12.4)') '  order ', order, '/', max_reflection_order, ' done, ratio = ', ratios(order)
      end if

      if (stop_ratio > 0.0_dp .and. ratios(order) < stop_ratio) then
        if (do_show) write(*,'(A,ES12.4,A)') '  stop_ratio reached: ratio = ', ratios(order), ' -> stop'
        exit
      end if
    end do

    deallocate(Jx_new, Jy_new, Jz_new)
  end subroutine build_iterative_po_currents_fast

  subroutine sum_currents_by_order(currents_x, currents_y, currents_z, n_orders, Jx_sum, Jy_sum, Jz_sum)
    complex(dp), intent(in) :: currents_x(:,:), currents_y(:,:), currents_z(:,:)
    integer, intent(in) :: n_orders
    complex(dp), intent(out) :: Jx_sum(:), Jy_sum(:), Jz_sum(:)

    Jx_sum = sum(currents_x(:,1:n_orders), dim=2)
    Jy_sum = sum(currents_y(:,1:n_orders), dim=2)
    Jz_sum = sum(currents_z(:,1:n_orders), dim=2)
  end subroutine sum_currents_by_order

end module mod_currents_solver
