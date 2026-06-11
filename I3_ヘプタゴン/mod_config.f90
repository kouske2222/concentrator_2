module mod_config
  use mod_types
  implicit none
  private
  public :: sim_config_type, load_default_config

  type :: sim_config_type
    ! ---- physical constants
    real(dp) :: c0, mu0, eps0, eta0

    ! ---- source / frequency
    real(dp) :: freq, lam0, k0
    real(dp) :: w0_src, z0_src, zR_src

    ! ---- geometry
    real(dp) :: d_in, base_z
    real(dp) :: r_cone_in, r_pipe
    real(dp) :: l_cone, l_pipe
    real(dp) :: f_cap
    real(dp) :: z_pipe_end, cap_depth, z_cap_vertex

    ! ---- mesh
    integer :: N_z_cone
    integer :: N_z_pipe
    integer :: N_r_cap
    integer :: N_theta_max

    ! ---- IPO
    integer :: max_reflection_order
    real(dp) :: stop_ratio
    real(dp) :: exclusion_factor_wall

    ! ---- field maps
    logical :: RUN_FIELD_MAPS
    logical :: RUN_XZ
    logical :: RUN_XY
    integer :: NX_XZ, NZ_XZ, BLOCK_SIZE_XZ
    integer :: NX_XY, NY_XY, N_XY_PLANES, BLOCK_SIZE_XY
    integer :: N_XY_PLANES_CONE
    integer :: N_XY_PLANES_PIPE 
    real(dp) :: exclusion_factor_obs
    real(dp) :: field_outer_margin_factor_xz
    real(dp) :: field_outer_margin_factor_xy

    ! ---- output
    character(len=256) :: OUTPUT_ROOT
    character(len=512) :: OUTPUT_DIR

    ! ---- progress / checkpoint
    logical :: SHOW_IPO_PROGRESS
    integer :: IPO_PROGRESS_UPDATES
    character(len=256) :: IPO_CHECKPOINT_DIR
    logical :: RESTART_FROM_CHECKPOINT
  end type sim_config_type

contains

  subroutine load_default_config(cfg)
    type(sim_config_type), intent(out) :: cfg

    cfg%c0   = 2.99792458e8_dp
    cfg%mu0  = 4.0e-7_dp * PI
    cfg%eps0 = 1.0_dp / (cfg%mu0 * cfg%c0 * cfg%c0)
    cfg%eta0 = sqrt(cfg%mu0 / cfg%eps0)

    cfg%freq = 94.0e9_dp
    cfg%lam0 = cfg%c0 / cfg%freq
    cfg%k0   = 2.0_dp * PI / cfg%lam0

    cfg%w0_src = 0.0204_dp
    cfg%z0_src = 0.0_dp
    cfg%zR_src = PI * cfg%w0_src**2 / cfg%lam0

    cfg%d_in      = 2.2_dp
    cfg%base_z    = cfg%d_in
    cfg%r_cone_in = 0.090_dp
    cfg%r_pipe    = 0.0280_dp
    cfg%l_cone    = 0.160_dp
    cfg%l_pipe    = 0.050_dp

    ! Parabolic closed cap:
    ! z(r) = z_pipe_end + (r_pipe^2 - r^2)/(4 f_cap)
    ! The cap opens toward -z and closes the +z end of the tube.
    cfg%f_cap        = 0.014_dp
    cfg%z_pipe_end   = cfg%base_z + cfg%l_cone + cfg%l_pipe
    cfg%cap_depth    = cfg%r_pipe**2 / (4.0_dp * cfg%f_cap)
    cfg%z_cap_vertex = cfg%z_pipe_end + cfg%cap_depth

    cfg%N_z_cone    = 600
    cfg%N_z_pipe    = 180
    cfg%N_r_cap     = 100
    cfg%N_theta_max = 1280

    cfg%max_reflection_order = 15
    cfg%stop_ratio            = 1.0e-2_dp
    cfg%exclusion_factor_wall = 0.15_dp

    cfg%RUN_FIELD_MAPS = .true.
    cfg%RUN_XZ = .true.
    cfg%RUN_XY = .true.

    cfg%NX_XZ = 400
    cfg%NZ_XZ = 500
    cfg%BLOCK_SIZE_XZ = 10000

    cfg%NX_XY = 401
    cfg%NY_XY = 401

    cfg%N_XY_PLANES_CONE = 8
    cfg%N_XY_PLANES_PIPE = 3
    cfg%N_XY_PLANES      = cfg%N_XY_PLANES_CONE + cfg%N_XY_PLANES_PIPE
    cfg%BLOCK_SIZE_XY = 10000

    cfg%exclusion_factor_obs = 0.15_dp
    cfg%field_outer_margin_factor_xz = 0.05_dp
    cfg%field_outer_margin_factor_xy = 0.05_dp

    cfg%OUTPUT_ROOT = 'results_heptagonal'
    cfg%OUTPUT_DIR  = ''

    cfg%SHOW_IPO_PROGRESS = .true.
    cfg%IPO_PROGRESS_UPDATES = 20

    ! The old open-end checkpoint is not compatible with the closed-cap geometry.
    cfg%IPO_CHECKPOINT_DIR = 'ipo_checkpoint_closed_cap_heptagonal'
    cfg%RESTART_FROM_CHECKPOINT = .true.
  end subroutine load_default_config

end module mod_config
