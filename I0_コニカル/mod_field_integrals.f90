module mod_field_integrals
  use mod_types
  implicit none
  private
  public :: calc_eh_from_surface_currents_to_points
  public :: calc_e_from_surface_currents_to_points

contains

  subroutine calc_eh_from_surface_currents_to_points( &
      src_Jx, src_Jy, src_Jz, src_x, src_y, src_z, src_dS, src_panel_len, &
      obs_x, obs_y, obs_z, k0, eta0, exclusion_factor, Ex, Ey, Ez, Hx, Hy, Hz, obs_mask)
    complex(dp), intent(in) :: src_Jx(:), src_Jy(:), src_Jz(:)
    real(dp), intent(in) :: src_x(:), src_y(:), src_z(:), src_dS(:), src_panel_len(:)
    real(dp), intent(in) :: obs_x(:), obs_y(:), obs_z(:)
    real(dp), intent(in) :: k0, eta0, exclusion_factor
    complex(dp), intent(out) :: Ex(:), Ey(:), Ez(:), Hx(:), Hy(:), Hz(:)
    logical, intent(in), optional :: obs_mask(:)

    integer :: i, s, N_src, N_obs
    real(dp) :: xo, yo, zo, rx, ry, rz, R, invR, ux, uy, uz
    real(dp) :: pref_e, pref_h
    complex(dp) :: ex_sum, ey_sum, ez_sum, hx_sum, hy_sum, hz_sum
    complex(dp) :: jx, jy, jz, ju, uju_x, uju_y, uju_z, phase, a, b, ch
    complex(dp) :: uxJx, uxJy, uxJz

    N_src = size(src_x)
    N_obs = size(obs_x)
    pref_e = eta0 / (4.0_dp * PI)
    pref_h = 1.0_dp / (4.0_dp * PI)

!$omp parallel do schedule(static) default(shared) &
!$omp private(i,s,xo,yo,zo,rx,ry,rz,R,invR,ux,uy,uz, &
!$omp ex_sum,ey_sum,ez_sum,hx_sum,hy_sum,hz_sum, &
!$omp jx,jy,jz,ju,uju_x,uju_y,uju_z,phase,a,b,ch,uxJx,uxJy,uxJz)
    do i = 1, N_obs
      if (present(obs_mask)) then
        if (.not. obs_mask(i)) then
          Ex(i) = (0.0_dp, 0.0_dp)
          Ey(i) = (0.0_dp, 0.0_dp)
          Ez(i) = (0.0_dp, 0.0_dp)
          Hx(i) = (0.0_dp, 0.0_dp)
          Hy(i) = (0.0_dp, 0.0_dp)
          Hz(i) = (0.0_dp, 0.0_dp)
          cycle
        end if
      end if

      xo = obs_x(i)
      yo = obs_y(i)
      zo = obs_z(i)

      ex_sum = (0.0_dp, 0.0_dp)
      ey_sum = (0.0_dp, 0.0_dp)
      ez_sum = (0.0_dp, 0.0_dp)
      hx_sum = (0.0_dp, 0.0_dp)
      hy_sum = (0.0_dp, 0.0_dp)
      hz_sum = (0.0_dp, 0.0_dp)

      do s = 1, N_src
        rx = xo - src_x(s)
        ry = yo - src_y(s)
        rz = zo - src_z(s)
        R = sqrt(rx*rx + ry*ry + rz*rz)
        if (R < exclusion_factor * src_panel_len(s)) cycle

        invR = 1.0_dp / R
        ux = rx * invR
        uy = ry * invR
        uz = rz * invR

        jx = src_Jx(s)
        jy = src_Jy(s)
        jz = src_Jz(s)

        ju = jx*ux + jy*uy + jz*uz
        uju_x = ux * ju
        uju_y = uy * ju
        uju_z = uz * ju

        phase = exp(-I_C*k0*R) * src_dS(s)
        a = (I_C*k0) * invR
        b = invR*invR - I_C*invR*invR*invR/k0

        ex_sum = ex_sum + phase * (a*(uju_x - jx) + b*(3.0_dp*uju_x - jx))
        ey_sum = ey_sum + phase * (a*(uju_y - jy) + b*(3.0_dp*uju_y - jy))
        ez_sum = ez_sum + phase * (a*(uju_z - jz) + b*(3.0_dp*uju_z - jz))

        uxJx = uy*jz - uz*jy
        uxJy = uz*jx - ux*jz
        uxJz = ux*jy - uy*jx
        ch = -(I_C*k0*invR + invR*invR)

        hx_sum = hx_sum + phase * ch * uxJx
        hy_sum = hy_sum + phase * ch * uxJy
        hz_sum = hz_sum + phase * ch * uxJz
      end do

      Ex(i) = pref_e * ex_sum
      Ey(i) = pref_e * ey_sum
      Ez(i) = pref_e * ez_sum
      Hx(i) = pref_h * hx_sum
      Hy(i) = pref_h * hy_sum
      Hz(i) = pref_h * hz_sum
    end do
!$omp end parallel do
  end subroutine calc_eh_from_surface_currents_to_points

  subroutine calc_e_from_surface_currents_to_points( &
      src_Jx, src_Jy, src_Jz, src_x, src_y, src_z, src_dS, src_panel_len, &
      obs_x, obs_y, obs_z, k0, eta0, exclusion_factor, Ex, Ey, Ez, obs_mask)
    complex(dp), intent(in) :: src_Jx(:), src_Jy(:), src_Jz(:)
    real(dp), intent(in) :: src_x(:), src_y(:), src_z(:), src_dS(:), src_panel_len(:)
    real(dp), intent(in) :: obs_x(:), obs_y(:), obs_z(:)
    real(dp), intent(in) :: k0, eta0, exclusion_factor
    complex(dp), intent(out) :: Ex(:), Ey(:), Ez(:)
    logical, intent(in), optional :: obs_mask(:)

    integer :: i, s, N_src, N_obs
    real(dp) :: xo, yo, zo, rx, ry, rz, R, invR, ux, uy, uz
    real(dp) :: pref_e
    complex(dp) :: ex_sum, ey_sum, ez_sum
    complex(dp) :: jx, jy, jz, ju, uju_x, uju_y, uju_z, phase, a, b

    N_src = size(src_x)
    N_obs = size(obs_x)
    pref_e = eta0 / (4.0_dp * PI)

!$omp parallel do schedule(static) default(shared) &
!$omp private(i,s,xo,yo,zo,rx,ry,rz,R,invR,ux,uy,uz, &
!$omp ex_sum,ey_sum,ez_sum,jx,jy,jz,ju,uju_x,uju_y,uju_z,phase,a,b)
    do i = 1, N_obs
      if (present(obs_mask)) then
        if (.not. obs_mask(i)) then
          Ex(i) = (0.0_dp, 0.0_dp)
          Ey(i) = (0.0_dp, 0.0_dp)
          Ez(i) = (0.0_dp, 0.0_dp)
          cycle
        end if
      end if

      xo = obs_x(i)
      yo = obs_y(i)
      zo = obs_z(i)

      ex_sum = (0.0_dp, 0.0_dp)
      ey_sum = (0.0_dp, 0.0_dp)
      ez_sum = (0.0_dp, 0.0_dp)

      do s = 1, N_src
        rx = xo - src_x(s)
        ry = yo - src_y(s)
        rz = zo - src_z(s)
        R = sqrt(rx*rx + ry*ry + rz*rz)
        if (R < exclusion_factor * src_panel_len(s)) cycle

        invR = 1.0_dp / R
        ux = rx * invR
        uy = ry * invR
        uz = rz * invR

        jx = src_Jx(s)
        jy = src_Jy(s)
        jz = src_Jz(s)

        ju = jx*ux + jy*uy + jz*uz
        uju_x = ux * ju
        uju_y = uy * ju
        uju_z = uz * ju

        phase = exp(-I_C*k0*R) * src_dS(s)
        a = (I_C*k0) * invR
        b = invR*invR - I_C*invR*invR*invR/k0

        ex_sum = ex_sum + phase * (a*(uju_x - jx) + b*(3.0_dp*uju_x - jx))
        ey_sum = ey_sum + phase * (a*(uju_y - jy) + b*(3.0_dp*uju_y - jy))
        ez_sum = ez_sum + phase * (a*(uju_z - jz) + b*(3.0_dp*uju_z - jz))
      end do

      Ex(i) = pref_e * ex_sum
      Ey(i) = pref_e * ey_sum
      Ez(i) = pref_e * ez_sum
    end do
!$omp end parallel do
  end subroutine calc_e_from_surface_currents_to_points

end module mod_field_integrals
