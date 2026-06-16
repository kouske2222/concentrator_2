module mod_geometry
  use mod_types
  implicit none
  private
  public :: wall_radius_at_z
  public :: wall_radius_at_theta_z
  public :: wall_radius_at_xy_z
  public :: cap_surface_z_of_r
  public :: octagon_softmax_beta_default
  public :: define_closed_concentrator_panels_adaptive

contains

  pure real(dp) function cap_surface_z_of_r(r, r_pipe, z_pipe_end, f_cap) result(zs)
    real(dp), intent(in) :: r, r_pipe, z_pipe_end, f_cap
    zs = z_pipe_end + (r_pipe*r_pipe - r*r) / (4.0_dp * f_cap)
  end function cap_surface_z_of_r

  pure real(dp) function wall_radius_at_z( &
      z, base_z, l_cone, l_pipe, r_cone_in, r_pipe, f_cap) result(rw)
    real(dp), intent(in) :: z, base_z, l_cone, l_pipe
    real(dp), intent(in) :: r_cone_in, r_pipe, f_cap
    real(dp) :: zloc, z_pipe_end, z_cap_vertex, arg
    real(dp) :: t, s

    zloc = z - base_z
    z_pipe_end = base_z + l_cone + l_pipe
    z_cap_vertex = z_pipe_end + r_pipe*r_pipe / (4.0_dp * f_cap)

    if (z < base_z .or. z > z_cap_vertex) then
      rw = -1.0_dp
    else if (zloc <= l_cone) then
      t = zloc / l_cone
      s = delayed_oct_to_circle_weight(t)

      ! Plot / grid-size radius for tapered octagon-to-circle geometry.
      rw = (1.0_dp - s) * octagonal_taper_radius(t, r_cone_in, r_pipe) + s * r_pipe
    else if (z <= z_pipe_end) then
      rw = r_pipe
    else
      arg = r_pipe*r_pipe - 4.0_dp * f_cap * (z - z_pipe_end)
      rw = sqrt(max(arg, 0.0_dp))
    end if
  end function wall_radius_at_z

  pure real(dp) function octagon_softmax_beta_default() result(beta)

    beta = 1.0_dp

  end function octagon_softmax_beta_default


  pure real(dp) function smootherstep_local(t) result(s)
    real(dp), intent(in) :: t
    real(dp) :: tt

    tt = min(max(t, 0.0_dp), 1.0_dp)

    s = 6.0_dp*tt**5 - 15.0_dp*tt**4 + 10.0_dp*tt**3

  end function smootherstep_local

  pure real(dp) function octagonal_taper_radius(t, r_cone_in, r_pipe) result(rref)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: r_cone_in
    real(dp), intent(in) :: r_pipe

    real(dp) :: tt

    tt = min(max(t, 0.0_dp), 1.0_dp)

    ! Linear taper in z.
    ! t = 0 : entrance octagonal radius
    ! t = 1 : pipe radius
    rref = r_cone_in + (r_pipe - r_cone_in) * tt

  end function octagonal_taper_radius

  ! Delayed transition from soft octagonal prism to circular pipe.
  ! t = 0.0 ... 1.0 corresponds to the cone section.
  ! For t <= oct_hold_fraction_default(), the cross-section remains a
  ! straight soft-octagonal prism.  Between t0 and t1 it changes smoothly
  ! to a circle.  At both ends, ds/dt = d2s/dt2 = 0.
  pure real(dp) function oct_hold_fraction_default() result(t0)

    t0 = 0.8_dp

  end function oct_hold_fraction_default


  pure real(dp) function oct_transition_end_fraction_default() result(t1)

    t1 = 0.95_dp

  end function oct_transition_end_fraction_default


  pure real(dp) function delayed_oct_to_circle_weight(t) result(s)
    real(dp), intent(in) :: t

    real(dp) :: t0, t1, u

    t0 = oct_hold_fraction_default()
    t1 = oct_transition_end_fraction_default()

    if (t <= t0) then
      s = 0.0_dp
    else if (t >= t1) then
      s = 1.0_dp
    else if (t1 <= t0) then
      s = 1.0_dp
    else
      u = (t - t0) / (t1 - t0)
      s = smootherstep_local(u)
    end if

  end function delayed_oct_to_circle_weight


  pure real(dp) function soft_octagon_radius(theta, r_vertex) result(r8)
    real(dp), intent(in) :: theta
    real(dp), intent(in) :: r_vertex

    integer :: k
    real(dp) :: beta
    real(dp) :: pi_l
    real(dp) :: alpha
    real(dp) :: cval
    real(dp) :: m_theta
    real(dp) :: m_vertex
    real(dp) :: sum_exp
    real(dp) :: max_arg
    real(dp) :: arg
    real(dp) :: theta0

    pi_l = acos(-1.0_dp)
    beta = octagon_softmax_beta_default()

    theta0 = modulo(theta, 2.0_dp*pi_l)

    ! ---- m(theta): softmax of cos(theta - alpha_k)
    max_arg = -huge(1.0_dp)
    do k = 0, 7
      alpha = pi_l/8.0_dp + real(k, dp)*pi_l/4.0_dp
      arg = beta*cos(theta0 - alpha)
      if (arg > max_arg) max_arg = arg
    end do

    sum_exp = 0.0_dp
    do k = 0, 7
      alpha = pi_l/8.0_dp + real(k, dp)*pi_l/4.0_dp
      arg = beta*cos(theta0 - alpha)
      sum_exp = sum_exp + exp(arg - max_arg)
    end do

    m_theta = (max_arg + log(sum_exp)) / beta

    ! ---- m(0): normalize so that radius at vertex direction is r_vertex
    max_arg = -huge(1.0_dp)
    do k = 0, 7
      alpha = pi_l/8.0_dp + real(k, dp)*pi_l/4.0_dp
      arg = beta*cos(0.0_dp - alpha)
      if (arg > max_arg) max_arg = arg
    end do

    sum_exp = 0.0_dp
    do k = 0, 7
      alpha = pi_l/8.0_dp + real(k, dp)*pi_l/4.0_dp
      arg = beta*cos(0.0_dp - alpha)
      sum_exp = sum_exp + exp(arg - max_arg)
    end do

    m_vertex = (max_arg + log(sum_exp)) / beta

    r8 = r_vertex * m_vertex / m_theta

  end function soft_octagon_radius


  pure real(dp) function blended_cone_radius(theta, t, r_cone_in, r_pipe) result(rw)
      real(dp), intent(in) :: theta
      real(dp), intent(in) :: t
      real(dp), intent(in) :: r_cone_in
      real(dp), intent(in) :: r_pipe

      real(dp) :: s
      real(dp) :: r8
      real(dp) :: rref

      ! Octagon-to-circle transition weight.
      ! s = 0 : octagonal cross-section
      ! s = 1 : circular cross-section
      s = delayed_oct_to_circle_weight(t)

      ! z-dependent reference radius.
      ! This makes the octagonal part a tapered frustum, not a prism.
      rref = octagonal_taper_radius(t, r_cone_in, r_pipe)

      ! Soft octagon whose size decreases along z.
      r8 = soft_octagon_radius(theta, rref)

      ! Blend from tapered octagon to circular pipe.
      rw = (1.0_dp - s) * r8 + s * r_pipe

  end function blended_cone_radius


  pure real(dp) function wall_radius_at_theta_z( &
      theta, z, base_z, l_cone, l_pipe, r_cone_in, r_pipe, f_cap) result(rw)

    real(dp), intent(in) :: theta
    real(dp), intent(in) :: z
    real(dp), intent(in) :: base_z, l_cone, l_pipe
    real(dp), intent(in) :: r_cone_in, r_pipe, f_cap

    real(dp) :: zloc
    real(dp) :: z_pipe_end
    real(dp) :: z_cap_vertex
    real(dp) :: arg
    real(dp) :: t
    real(dp) :: theta_use

    zloc = z - base_z
    z_pipe_end = base_z + l_cone + l_pipe
    z_cap_vertex = z_pipe_end + r_pipe*r_pipe / (4.0_dp * f_cap)

    theta_use = modulo(theta, 2.0_dp * PI)

    if (z < base_z .or. z > z_cap_vertex) then
      rw = -1.0_dp

    else if (zloc <= l_cone) then
      t = zloc / l_cone
      rw = blended_cone_radius(theta_use, t, r_cone_in, r_pipe)

    else if (z <= z_pipe_end) then
      rw = r_pipe

    else
      arg = r_pipe*r_pipe - 4.0_dp * f_cap * (z - z_pipe_end)
      rw = sqrt(max(arg, 0.0_dp))
    end if

  end function wall_radius_at_theta_z


  pure real(dp) function wall_radius_at_xy_z( &
      x, y, z, base_z, l_cone, l_pipe, r_cone_in, r_pipe, f_cap) result(rw)

    real(dp), intent(in) :: x, y, z
    real(dp), intent(in) :: base_z, l_cone, l_pipe
    real(dp), intent(in) :: r_cone_in, r_pipe, f_cap

    real(dp) :: theta

    if (x == 0.0_dp .and. y == 0.0_dp) then
      theta = 0.0_dp
    else
      theta = atan2(y, x)
      if (theta < 0.0_dp) theta = theta + 2.0_dp * PI
    end if

    rw = wall_radius_at_theta_z( &
        theta, z, base_z, l_cone, l_pipe, r_cone_in, r_pipe, f_cap)

  end function wall_radius_at_xy_z

  subroutine define_closed_concentrator_panels_adaptive( &
      r_cone_in, r_pipe, l_cone, l_pipe, base_z, f_cap, &
      N_z_cone, N_z_pipe, N_r_cap, N_theta_max, &
      X, Y, Z, nx, ny, nz, dS, zone_ids)
    real(dp), intent(in) :: r_cone_in, r_pipe, l_cone, l_pipe, base_z, f_cap
    integer, intent(in) :: N_z_cone, N_z_pipe, N_r_cap, N_theta_max
    real(dp), allocatable, intent(out) :: X(:), Y(:), Z(:), nx(:), ny(:), nz(:), dS(:)
    integer, allocatable, intent(out) :: zone_ids(:)

    integer :: i, j, idx, n_total
    integer :: N_th, N_th_p
    real(dp) :: ds_target
    real(dp) :: dz_c, dz_p, dr_cap
    real(dp) :: norm_cap
    real(dp) :: zc, zp, Rp
    real(dp) :: r_mid, z_cap, dtheta, dtheta_p, th, cth, sth
    real(dp) :: z_pipe_end
    real(dp) :: t_c, s_c, r_count
    real(dp) :: r_here, r_th_p, r_th_m, r_z_p, r_z_m
    real(dp) :: eps_th, eps_z, dr_dth, dr_dzc
    real(dp) :: eth_x, eth_y, eth_z
    real(dp) :: ez_x, ez_y, ez_z
    real(dp) :: cx, cy, cz, cnorm
    integer, allocatable :: ntheta_cone(:), ntheta_cap(:)

    ds_target = 2.0_dp * PI * r_cone_in / real(N_theta_max, dp)

    ! ---- count cone panels
    dz_c = l_cone / real(N_z_cone, dp)

    allocate(ntheta_cone(N_z_cone))
    n_total = 0
    do i = 1, N_z_cone
      zc = (real(i, dp) - 0.5_dp) * dz_c
      t_c = zc / l_cone
      s_c = delayed_oct_to_circle_weight(t_c)
      ! Maximum radius at this z.  The actual radius still depends on theta.
      r_count = (1.0_dp - s_c) * octagonal_taper_radius(t_c, r_cone_in, r_pipe)  + s_c * r_pipe
      ntheta_cone(i) = max(nint(2.0_dp * PI * r_count / ds_target), 12)
      n_total = n_total + ntheta_cone(i)
    end do

    ! ---- count pipe panels
    dz_p = l_pipe / real(N_z_pipe, dp)
    Rp = r_pipe
    N_th_p = max(nint(2.0_dp * PI * Rp / ds_target), 12)
    n_total = n_total + N_z_pipe * N_th_p

    ! ---- count paraboloid-cap panels
    dr_cap = r_pipe / real(N_r_cap, dp)
    allocate(ntheta_cap(N_r_cap))
    do i = 1, N_r_cap
      r_mid = (real(i, dp) - 0.5_dp) * dr_cap
      ntheta_cap(i) = max(nint(2.0_dp * PI * r_mid / ds_target), 12)
      n_total = n_total + ntheta_cap(i)
    end do

    allocate(X(n_total), Y(n_total), Z(n_total))
    allocate(nx(n_total), ny(n_total), nz(n_total), dS(n_total))
    allocate(zone_ids(n_total))

    idx = 0

    ! ============================================================
    ! zone 0: soft-octagonal prism followed by G2 transition to circular pipe
    !
    ! Surface parameterization:
    !   S(theta,z) = [r(theta,t) cos(theta),
    !                 r(theta,t) sin(theta),
    !                 base_z + z]
    !   t = z / l_cone; transition starts after oct_hold_fraction_default()
    ! The inward normal is -normalize(dS/dtheta x dS/dz).
    ! ============================================================
    eps_z = min(1.0e-5_dp * l_cone, 0.25_dp * dz_c)

    do i = 1, N_z_cone
      zc = (real(i, dp) - 0.5_dp) * dz_c
      t_c = zc / l_cone
      N_th = ntheta_cone(i)
      dtheta = 2.0_dp * PI / real(N_th, dp)
      eps_th = min(1.0e-5_dp, 0.25_dp * dtheta)

      do j = 1, N_th
        th = (real(j, dp) - 0.5_dp) * dtheta
        cth = cos(th)
        sth = sin(th)

        r_here = blended_cone_radius(th, t_c, r_cone_in, r_pipe)

        ! Numerical derivatives of r(theta,z).
        r_th_p = blended_cone_radius(th + eps_th, t_c, r_cone_in, r_pipe)
        r_th_m = blended_cone_radius(th - eps_th, t_c, r_cone_in, r_pipe)
        dr_dth = (r_th_p - r_th_m) / (2.0_dp * eps_th)

        r_z_p = blended_cone_radius(th, (zc + eps_z) / l_cone, r_cone_in, r_pipe)
        r_z_m = blended_cone_radius(th, (zc - eps_z) / l_cone, r_cone_in, r_pipe)
        dr_dzc = (r_z_p - r_z_m) / (2.0_dp * eps_z)

        ! Tangent vectors dS/dtheta and dS/dz.
        eth_x = dr_dth * cth - r_here * sth
        eth_y = dr_dth * sth + r_here * cth
        eth_z = 0.0_dp

        ez_x = dr_dzc * cth
        ez_y = dr_dzc * sth
        ez_z = 1.0_dp

        ! outward area vector = dS/dtheta x dS/dz
        cx = eth_y * ez_z - eth_z * ez_y
        cy = eth_z * ez_x - eth_x * ez_z
        cz = eth_x * ez_y - eth_y * ez_x
        cnorm = sqrt(cx*cx + cy*cy + cz*cz)

        idx = idx + 1
        X(idx) = r_here * cth
        Y(idx) = r_here * sth
        Z(idx) = base_z + zc

        ! inward normal, i.e. toward the cavity
        nx(idx) = -cx / cnorm
        ny(idx) = -cy / cnorm
        nz(idx) = -cz / cnorm

        dS(idx) = cnorm * dtheta * dz_c
        zone_ids(idx) = 0
      end do
    end do

    ! ============================================================
    ! zone 1: straight pipe
    ! ============================================================
    do i = 1, N_z_pipe
      zp = (real(i, dp) - 0.5_dp) * dz_p
      dtheta_p = 2.0_dp * PI / real(N_th_p, dp)

      do j = 1, N_th_p
        th = (real(j, dp) - 0.5_dp) * dtheta_p
        cth = cos(th)
        sth = sin(th)

        idx = idx + 1
        X(idx) = Rp * cth
        Y(idx) = Rp * sth
        Z(idx) = base_z + l_cone + zp

        nx(idx) = -cth
        ny(idx) = -sth
        nz(idx) =  0.0_dp

        dS(idx) = Rp * dtheta_p * dz_p
        zone_ids(idx) = 1
      end do
    end do

    ! ============================================================
    ! zone 2: rotational paraboloid closing the +z end
    !
    ! z(r) = z_pipe_end + (r_pipe^2 - r^2)/(4 f_cap)
    ! inward normal = -grad(F)/|grad(F)|, F=z-z(r)
    ! ============================================================
    z_pipe_end = base_z + l_cone + l_pipe

    do i = 1, N_r_cap
      r_mid = (real(i, dp) - 0.5_dp) * dr_cap
      z_cap = cap_surface_z_of_r(r_mid, r_pipe, z_pipe_end, f_cap)

      N_th = ntheta_cap(i)
      dtheta = 2.0_dp * PI / real(N_th, dp)
      norm_cap = sqrt(1.0_dp + (r_mid / (2.0_dp * f_cap))**2)

      do j = 1, N_th
        th = (real(j, dp) - 0.5_dp) * dtheta
        cth = cos(th)
        sth = sin(th)

        idx = idx + 1
        X(idx) = r_mid * cth
        Y(idx) = r_mid * sth
        Z(idx) = z_cap

        nx(idx) = -(r_mid * cth / (2.0_dp * f_cap)) / norm_cap
        ny(idx) = -(r_mid * sth / (2.0_dp * f_cap)) / norm_cap
        nz(idx) = -1.0_dp / norm_cap

        dS(idx) = r_mid * norm_cap * dtheta * dr_cap
        zone_ids(idx) = 2
      end do
    end do

    if (idx /= n_total) stop 'Panel count mismatch in define_closed_concentrator_panels_adaptive.'

    deallocate(ntheta_cone, ntheta_cap)
  end subroutine define_closed_concentrator_panels_adaptive

end module mod_geometry
