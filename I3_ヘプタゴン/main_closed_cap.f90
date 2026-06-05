program main_closed_cap
  use mod_types
  use mod_config
  use mod_geometry
  use mod_incident
  use mod_field_integrals
  use mod_currents_solver
  use mod_io
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none

  type(sim_config_type) :: cfg

  real(dp), allocatable :: Xm(:), Ym(:), Zm(:)
  real(dp), allocatable :: nxm(:), nym(:), nzm(:)
  real(dp), allocatable :: dSm(:), panel_len_m(:)
  integer, allocatable :: zone_ids(:)

  complex(dp), allocatable :: currents_x(:,:), currents_y(:,:), currents_z(:,:)
  complex(dp), allocatable :: Jx_tot(:), Jy_tot(:), Jz_tot(:)
  real(dp), allocatable :: ratios(:), order_build_times(:)

  real(dp), allocatable :: obs_x(:), obs_y(:), obs_z(:), e_sq(:)
  complex(dp), allocatable :: Ex(:), Ey(:), Ez(:)
  complex(dp), allocatable :: Ex_s(:), Ey_s(:), Ez_s(:)
  logical, allocatable :: obs_mask(:)

  real(dp), allocatable :: xy_plane_z_list(:)

  integer :: Np, n_orders_done
  integer :: n_cone, n_pipe, n_cap
  integer :: Nx, Ny, Nz, ix, iy, iz, izp, idx
  integer :: n_total_xz, n_total_xy, ibeg, iend
  real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max
  real(dp) :: dx, dy, dz
  real(dp) :: z_plane, r_wall_here, r_plot_here, x_half_width_xz
  real(dp) :: t0, t1
  character(len=128) :: fname_xy

  call load_default_config(cfg)
  call create_timestamp_output_dir(cfg%OUTPUT_ROOT, cfg%OUTPUT_DIR)
  call cpu_time(t0)

  print *, '========================================'
  print *, 'Microwave PO / IPO: closed parabolic cap, heptagon cone'
  print *, 'PTD: OFF'
  print *, 'Output directory: ', trim(cfg%OUTPUT_DIR)
  print *, '========================================'

  ! ============================================================
  ! Step 1: geometry
  ! ============================================================
  print *, 'Step 1: Generate heptagon cone + pipe + parabolic-cap panels'

  call define_closed_concentrator_panels_adaptive( &
      cfg%r_cone_in, cfg%r_pipe, cfg%l_cone, cfg%l_pipe, cfg%base_z, cfg%f_cap, &
      cfg%N_z_cone, cfg%N_z_pipe, cfg%N_r_cap, cfg%N_theta_max, &
      Xm, Ym, Zm, nxm, nym, nzm, dSm, zone_ids)

  Np = size(Xm)
  allocate(panel_len_m(Np))
  allocate(Jx_tot(Np), Jy_tot(Np), Jz_tot(Np))
  panel_len_m = sqrt(dSm / PI)

  n_cone = count(zone_ids == 0)
  n_pipe = count(zone_ids == 1)
  n_cap  = count(zone_ids == 2)

  print *, '  done'
  print *, '  Total panels          = ', Np
  print *, '  cone panels           = ', n_cone
  print *, '  pipe panels           = ', n_pipe
  print *, '  parabolic-cap panels  = ', n_cap
  print *, '  cap depth [m]         = ', cfg%cap_depth
  print *, '  cap vertex z [m]      = ', cfg%z_cap_vertex

  call export_geometry_summary( &
      trim(cfg%OUTPUT_DIR)//'/geometry_summary.txt', cfg%r_cone_in, cfg%r_pipe, cfg%l_cone, cfg%l_pipe, &
      cfg%base_z, cfg%f_cap, cfg%cap_depth, cfg%z_cap_vertex, &
      Np, n_cone, n_pipe, n_cap)

  ! ============================================================
  ! Step 2: IPO currents
  ! ============================================================
  print *, 'Step 2: Build iterative PO currents'

  call build_iterative_po_currents_fast( &
      Xm, Ym, Zm, nxm, nym, nzm, dSm, panel_len_m, &
      cfg%max_reflection_order, cfg%k0, cfg%eta0, cfg%w0_src, cfg%z0_src, cfg%zR_src, &
      cfg%exclusion_factor_wall, cfg%stop_ratio, &
      currents_x, currents_y, currents_z, ratios, order_build_times, n_orders_done, &
      show_progress=cfg%SHOW_IPO_PROGRESS, n_progress_updates=cfg%IPO_PROGRESS_UPDATES, &
      checkpoint_dir=trim(cfg%IPO_CHECKPOINT_DIR), restart_from_checkpoint=cfg%RESTART_FROM_CHECKPOINT)

  print *, '  done'
  print *, '  IPO orders used = ', n_orders_done

  call sum_currents_by_order(currents_x, currents_y, currents_z, n_orders_done, Jx_tot, Jy_tot, Jz_tot)
  call export_ipo_summary(trim(cfg%OUTPUT_DIR)//'/ipo_order_summary.txt', ratios, order_build_times, n_orders_done)

  ! ============================================================
  ! Step 3: field maps
  ! ============================================================
  if (cfg%RUN_FIELD_MAPS) then
    cfg%N_XY_PLANES = cfg%N_XY_PLANES_CONE + cfg%N_XY_PLANES_PIPE

    allocate(xy_plane_z_list(cfg%N_XY_PLANES))
    call build_split_xy_plane_list(cfg, xy_plane_z_list)
    call export_xy_plane_index_split( &
    trim(cfg%OUTPUT_DIR)//'/xy_plane_index_to_z.txt', &
    xy_plane_z_list, cfg%N_XY_PLANES_CONE)

    if (cfg%RUN_XZ) then
      call export_xz_field_map(cfg, Xm, Ym, Zm, dSm, panel_len_m, Jx_tot, Jy_tot, Jz_tot)
    end if

    if (cfg%RUN_XY) then
      call export_xy_field_maps(cfg, Xm, Ym, Zm, dSm, panel_len_m, Jx_tot, Jy_tot, Jz_tot, xy_plane_z_list)
    end if

    deallocate(xy_plane_z_list)
  else
    print *, 'Step 3: field maps skipped'
  end if

  call cpu_time(t1)

  print *, '========================================'
  print *, 'FINISHED'
  print *, 'Elapsed CPU time [s] = ', t1 - t0
  print *, 'Output directory:'
  print *, '  ', trim(cfg%OUTPUT_DIR)
  print *, 'Output files:'
  print *, '  geometry_summary.txt'
  print *, '  ipo_order_summary.txt'
  if (cfg%RUN_FIELD_MAPS) then
    if (cfg%RUN_XZ) print *, '  field_map_xz.bin'
    if (cfg%RUN_XY) then
      print *, '  field_map_xy_z###.bin'
      print *, '  xy_plane_index_to_z.txt'
    end if
  end if
  print *, '========================================'

contains

    subroutine build_split_xy_plane_list(cfg_in, z_list)
      type(sim_config_type), intent(in) :: cfg_in
      real(dp), intent(out) :: z_list(:)

      integer :: n_cone_l, n_pipe_l, n_total_l
      integer :: ip_l

      n_cone_l = max(cfg_in%N_XY_PLANES_CONE, 0)
      n_pipe_l = max(cfg_in%N_XY_PLANES_PIPE, 0)
      n_total_l = n_cone_l + n_pipe_l

      if (n_total_l <= 0) then
        stop 'N_XY_PLANES_CONE + N_XY_PLANES_PIPE must be positive.'
      end if

      if (size(z_list) /= n_total_l) then
        stop 'size(z_list) must equal N_XY_PLANES_CONE + N_XY_PLANES_PIPE.'
      end if

      ! ============================================================
      ! cone planes
      !   n_cone_l >= 2:
      !     include both cone inlet and cone outlet
      !   n_cone_l == 1:
      !     place at cone center
      ! ============================================================
      if (n_cone_l == 1) then
        z_list(1) = cfg_in%base_z + 0.5_dp * cfg_in%l_cone

      else if (n_cone_l >= 2) then
        do ip_l = 1, n_cone_l
          z_list(ip_l) = cfg_in%base_z &
              + cfg_in%l_cone * real(ip_l - 1, dp) / real(n_cone_l - 1, dp)
        end do
      end if

      ! ============================================================
      ! pipe planes
      !   avoid duplicating z = base_z + l_cone
      !   include pipe end
      ! ============================================================
      if (n_pipe_l == 1) then
        z_list(n_cone_l + 1) = cfg_in%base_z + cfg_in%l_cone &
            + 0.5_dp * cfg_in%l_pipe

      else if (n_pipe_l >= 2) then
        do ip_l = 1, n_pipe_l
          z_list(n_cone_l + ip_l) = cfg_in%base_z + cfg_in%l_cone &
              + cfg_in%l_pipe * real(ip_l, dp) / real(n_pipe_l, dp)
        end do
      end if

    end subroutine build_split_xy_plane_list

    subroutine export_xy_plane_index_split(filename, z_list, n_cone)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: z_list(:)
    integer, intent(in) :: n_cone

    integer :: u, i
    character(len=16) :: region

    open(newunit=u, file=trim(filename), status='replace', action='write')

    write(u,'(A)') '# index   z[m]                 region'

    do i = 1, size(z_list)
      if (i <= n_cone) then
        region = 'cone'
      else
        region = 'pipe'
      end if

      write(u,'(I6,1X,ES20.12,1X,A)') i, z_list(i), trim(region)
    end do

    close(u)
  end subroutine export_xy_plane_index_split


  subroutine export_xz_field_map(cfg_in, Xm_in, Ym_in, Zm_in, dSm_in, panel_len_in, Jx_in, Jy_in, Jz_in)
    type(sim_config_type), intent(in) :: cfg_in
    real(dp), intent(in) :: Xm_in(:), Ym_in(:), Zm_in(:), dSm_in(:), panel_len_in(:)
    complex(dp), intent(in) :: Jx_in(:), Jy_in(:), Jz_in(:)

    real(dp), allocatable :: xg(:), yg(:), zg(:), esq(:)
    complex(dp), allocatable :: Ex_tot(:), Ey_tot(:), Ez_tot(:)
    complex(dp), allocatable :: Ex_sca(:), Ey_sca(:), Ez_sca(:)
    logical, allocatable :: mask(:)

    integer :: Nx_l, Nz_l, ntot_l, ix_l, iz_l, idx_l, ibeg_l, iend_l
    real(dp) :: x_min_l, x_max_l, z_min_l, z_max_l, dx_l, dz_l
    real(dp) :: r_wall_l, x_half_l

    print *, 'Step 3a: Calculating XZ field map'

    Nx_l = cfg_in%NX_XZ
    Nz_l = cfg_in%NZ_XZ
    ntot_l = Nx_l * Nz_l

    allocate(xg(ntot_l), yg(ntot_l), zg(ntot_l), esq(ntot_l))
    allocate(mask(ntot_l))
    allocate(Ex_tot(ntot_l), Ey_tot(ntot_l), Ez_tot(ntot_l))
    allocate(Ex_sca(ntot_l), Ey_sca(ntot_l), Ez_sca(ntot_l))

    x_half_l = (1.0_dp + cfg_in%field_outer_margin_factor_xz) * max(cfg_in%r_cone_in, cfg_in%r_pipe)
    x_min_l = -x_half_l
    x_max_l =  x_half_l
    z_min_l = cfg_in%base_z
    z_max_l = cfg_in%z_cap_vertex

    dx_l = (x_max_l - x_min_l) / real(Nx_l - 1, dp)
    dz_l = (z_max_l - z_min_l) / real(Nz_l - 1, dp)

    idx_l = 0
    do iz_l = 1, Nz_l
      do ix_l = 1, Nx_l
        idx_l = idx_l + 1
        xg(idx_l) = x_min_l + real(ix_l - 1, dp) * dx_l
        yg(idx_l) = 0.0_dp
        zg(idx_l) = z_min_l + real(iz_l - 1, dp) * dz_l

        r_wall_l = wall_radius_at_z( &
            zg(idx_l), cfg_in%base_z, cfg_in%l_cone, cfg_in%l_pipe, &
            cfg_in%r_cone_in, cfg_in%r_pipe, cfg_in%f_cap)

        mask(idx_l) = (r_wall_l >= 0.0_dp) .and. &
                      (abs(xg(idx_l)) <= (1.0_dp + cfg_in%field_outer_margin_factor_xz) * r_wall_l)
      end do
    end do

    call calc_gaussian_incident_e_points( &
        xg, yg, zg, cfg_in%w0_src, cfg_in%z0_src, cfg_in%k0, cfg_in%zR_src, &
        Ex_tot, Ey_tot, Ez_tot)

    where (.not. mask)
      Ex_tot = (0.0_dp, 0.0_dp)
      Ey_tot = (0.0_dp, 0.0_dp)
      Ez_tot = (0.0_dp, 0.0_dp)
    end where

    Ex_sca = (0.0_dp, 0.0_dp)
    Ey_sca = (0.0_dp, 0.0_dp)
    Ez_sca = (0.0_dp, 0.0_dp)

    do ibeg_l = 1, ntot_l, cfg_in%BLOCK_SIZE_XZ
      iend_l = min(ibeg_l + cfg_in%BLOCK_SIZE_XZ - 1, ntot_l)
      call calc_e_from_surface_currents_to_points( &
          Jx_in, Jy_in, Jz_in, Xm_in, Ym_in, Zm_in, dSm_in, panel_len_in, &
          xg(ibeg_l:iend_l), yg(ibeg_l:iend_l), zg(ibeg_l:iend_l), &
          cfg_in%k0, cfg_in%eta0, cfg_in%exclusion_factor_obs, &
          Ex_sca(ibeg_l:iend_l), Ey_sca(ibeg_l:iend_l), Ez_sca(ibeg_l:iend_l), &
          mask(ibeg_l:iend_l))
    end do

    Ex_tot = Ex_tot + Ex_sca
    Ey_tot = Ey_tot + Ey_sca
    Ez_tot = Ez_tot + Ez_sca

    esq = ieee_value(0.0_dp, ieee_quiet_nan)
!$omp parallel do schedule(static) default(shared) private(idx_l)
    do idx_l = 1, ntot_l
      if (mask(idx_l)) then
        esq(idx_l) = abs(Ex_tot(idx_l))**2 + abs(Ey_tot(idx_l))**2 + abs(Ez_tot(idx_l))**2
      end if
    end do
!$omp end parallel do

    call export_field_map_binary(trim(cfg_in%OUTPUT_DIR)//'/field_map_xz.bin', ntot_l, xg, yg, zg, esq)
    print *, '  XZ field map exported: ', trim(cfg_in%OUTPUT_DIR)//'/field_map_xz.bin'

    deallocate(xg, yg, zg, esq, mask)
    deallocate(Ex_tot, Ey_tot, Ez_tot, Ex_sca, Ey_sca, Ez_sca)
  end subroutine export_xz_field_map

  subroutine export_xy_field_maps( &
      cfg_in, Xm_in, Ym_in, Zm_in, dSm_in, panel_len_in, Jx_in, Jy_in, Jz_in, z_list)
    type(sim_config_type), intent(in) :: cfg_in
    real(dp), intent(in) :: Xm_in(:), Ym_in(:), Zm_in(:), dSm_in(:), panel_len_in(:)
    complex(dp), intent(in) :: Jx_in(:), Jy_in(:), Jz_in(:)
    real(dp), intent(in) :: z_list(:)
    real(dp) :: z_l, r_wall_l, r_plot_l
real(dp) :: r_cell_l, theta_l

    real(dp), allocatable :: xg(:), yg(:), zg(:), esq(:)
    complex(dp), allocatable :: Ex_tot(:), Ey_tot(:), Ez_tot(:)
    complex(dp), allocatable :: Ex_sca(:), Ey_sca(:), Ez_sca(:)
    logical, allocatable :: mask(:)

    integer :: Nx_l, Ny_l, ntot_l
    integer :: ix_l, iy_l, izp_l, idx_l, ibeg_l, iend_l
    real(dp) :: x_min_l, x_max_l, y_min_l, y_max_l, dx_l, dy_l
    character(len=512) :: fname_l

    print *, 'Step 3b: Calculating XY field maps'

    Nx_l = cfg_in%NX_XY
    Ny_l = cfg_in%NY_XY
    ntot_l = Nx_l * Ny_l

    allocate(xg(ntot_l), yg(ntot_l), zg(ntot_l), esq(ntot_l))
    allocate(mask(ntot_l))
    allocate(Ex_tot(ntot_l), Ey_tot(ntot_l), Ez_tot(ntot_l))
    allocate(Ex_sca(ntot_l), Ey_sca(ntot_l), Ez_sca(ntot_l))

    do izp_l = 1, size(z_list)
      z_l = z_list(izp_l)

      r_wall_l = wall_radius_at_z( &
          z_l, cfg_in%base_z, cfg_in%l_cone, cfg_in%l_pipe, &
          cfg_in%r_cone_in, cfg_in%r_pipe, cfg_in%f_cap)

      if (r_wall_l < 0.0_dp) then
        write(*,'(A,I0,A)') '  skip XY plane ', izp_l, ': outside geometry'
        cycle
      end if

      r_plot_l = (1.0_dp + cfg_in%field_outer_margin_factor_xy) * r_wall_l
      x_min_l = -r_plot_l
      x_max_l =  r_plot_l
      y_min_l = -r_plot_l
      y_max_l =  r_plot_l

      dx_l = (x_max_l - x_min_l) / real(Nx_l - 1, dp)
      dy_l = (y_max_l - y_min_l) / real(Ny_l - 1, dp)

      idx_l = 0
        do iy_l = 1, Ny_l
          do ix_l = 1, Nx_l
            idx_l = idx_l + 1

            xg(idx_l) = x_min_l + real(ix_l - 1, dp) * dx_l
            yg(idx_l) = y_min_l + real(iy_l - 1, dp) * dy_l
            zg(idx_l) = z_l

            r_cell_l = sqrt(xg(idx_l)**2 + yg(idx_l)**2)

            r_wall_l = wall_radius_at_xy_z( &
                xg(idx_l), yg(idx_l), z_l, &
                cfg_in%base_z, cfg_in%l_cone, cfg_in%l_pipe, &
                cfg_in%r_cone_in, cfg_in%r_pipe, cfg_in%f_cap)

            mask(idx_l) = (r_wall_l >= 0.0_dp) .and. &
                          (r_cell_l <= (1.0_dp + cfg_in%field_outer_margin_factor_xy) * r_wall_l)
          end do
        end do

      call calc_gaussian_incident_e_points( &
          xg, yg, zg, cfg_in%w0_src, cfg_in%z0_src, cfg_in%k0, cfg_in%zR_src, &
          Ex_tot, Ey_tot, Ez_tot)

      where (.not. mask)
        Ex_tot = (0.0_dp, 0.0_dp)
        Ey_tot = (0.0_dp, 0.0_dp)
        Ez_tot = (0.0_dp, 0.0_dp)
      end where

      Ex_sca = (0.0_dp, 0.0_dp)
      Ey_sca = (0.0_dp, 0.0_dp)
      Ez_sca = (0.0_dp, 0.0_dp)

      do ibeg_l = 1, ntot_l, cfg_in%BLOCK_SIZE_XY
        iend_l = min(ibeg_l + cfg_in%BLOCK_SIZE_XY - 1, ntot_l)
        call calc_e_from_surface_currents_to_points( &
            Jx_in, Jy_in, Jz_in, Xm_in, Ym_in, Zm_in, dSm_in, panel_len_in, &
            xg(ibeg_l:iend_l), yg(ibeg_l:iend_l), zg(ibeg_l:iend_l), &
            cfg_in%k0, cfg_in%eta0, cfg_in%exclusion_factor_obs, &
            Ex_sca(ibeg_l:iend_l), Ey_sca(ibeg_l:iend_l), Ez_sca(ibeg_l:iend_l), &
            mask(ibeg_l:iend_l))
      end do

      Ex_tot = Ex_tot + Ex_sca
      Ey_tot = Ey_tot + Ey_sca
      Ez_tot = Ez_tot + Ez_sca

      esq = ieee_value(0.0_dp, ieee_quiet_nan)
!$omp parallel do schedule(static) default(shared) private(idx_l)
      do idx_l = 1, ntot_l
        if (mask(idx_l)) then
          esq(idx_l) = abs(Ex_tot(idx_l))**2 + abs(Ey_tot(idx_l))**2 + abs(Ez_tot(idx_l))**2
        end if
      end do
!$omp end parallel do

      write(fname_l,'(A,"/field_map_xy_z",I3.3,".bin")') trim(cfg_in%OUTPUT_DIR), izp_l
      call export_field_map_binary(trim(fname_l), ntot_l, xg, yg, zg, esq)
      write(*,'(A,I0,A,A)') '  XY plane ', izp_l, ' exported: ', trim(fname_l)
    end do

    deallocate(xg, yg, zg, esq, mask)
    deallocate(Ex_tot, Ey_tot, Ez_tot, Ex_sca, Ey_sca, Ez_sca)
  end subroutine export_xy_field_maps

end program main_closed_cap
