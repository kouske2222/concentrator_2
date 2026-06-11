module mod_incident
  use mod_types
  implicit none
  private
  public :: gaussian_scalar_field, gaussian_khat_and_polarization
  public :: calc_gaussian_incident_e_points
  public :: calc_gaussian_incident_eh_points
  public :: calc_surface_currents_from_incident_eh

contains

  pure complex(dp) function gaussian_scalar_field(x, y, z, w0, z0, k0, zR) result(E0)
    real(dp), intent(in) :: x, y, z, w0, z0, k0, zR
    real(dp) :: dz, r2, wz, psi, inv_R, Rz, amp, phase

    dz = z - z0
    r2 = x*x + y*y
    wz = w0 * sqrt(1.0_dp + (dz / zR)**2)
    psi = atan(dz / zR)

    if (abs(dz) < 1.0e-14_dp) then
      inv_R = 0.0_dp
    else
      Rz = dz * (1.0_dp + (zR / dz)**2)
      inv_R = 1.0_dp / Rz
    end if

    amp = (w0 / wz) * exp(-r2 / (wz*wz))
    phase = k0*dz + 0.5_dp*k0*r2*inv_R - psi
    E0 = amp * exp(-I_C * phase)
  end function gaussian_scalar_field

  pure subroutine gaussian_khat_and_polarization( &
      x, y, z, z0, zR, khx, khy, khz, px, py, pz)
    real(dp), intent(in) :: x, y, z, z0, zR
    real(dp), intent(out) :: khx, khy, khz, px, py, pz
    real(dp) :: dz, inv_R, Rz, kx, ky, kz, kn
    real(dp) :: ex0, ey0, ez0, dot0, ptx, pty, ptz, pn

    dz = z - z0
    if (abs(dz) < 1.0e-14_dp) then
      inv_R = 0.0_dp
    else
      Rz = dz * (1.0_dp + (zR / dz)**2)
      inv_R = 1.0_dp / Rz
    end if

    kx = x * inv_R
    ky = y * inv_R
    kz = 1.0_dp
    kn = sqrt(kx*kx + ky*ky + kz*kz)
    khx = kx / kn
    khy = ky / kn
    khz = kz / kn

    ex0 = 1.0_dp
    ey0 = 0.0_dp
    ez0 = 0.0_dp
    dot0 = ex0*khx + ey0*khy + ez0*khz

    ptx = ex0 - dot0*khx
    pty = ey0 - dot0*khy
    ptz = ez0 - dot0*khz
    pn = sqrt(ptx*ptx + pty*pty + ptz*ptz)

    if (pn < 1.0e-12_dp) then
      ex0 = 0.0_dp
      ey0 = 1.0_dp
      ez0 = 0.0_dp
      dot0 = ex0*khx + ey0*khy + ez0*khz
      ptx = ex0 - dot0*khx
      pty = ey0 - dot0*khy
      ptz = ez0 - dot0*khz
      pn = sqrt(ptx*ptx + pty*pty + ptz*ptz)
    end if

    px = ptx / pn
    py = pty / pn
    pz = ptz / pn
  end subroutine gaussian_khat_and_polarization

  subroutine calc_gaussian_incident_e_points(obs_x, obs_y, obs_z, w0, z0, k0, zR, Ex, Ey, Ez)
    real(dp), intent(in) :: obs_x(:), obs_y(:), obs_z(:)
    real(dp), intent(in) :: w0, z0, k0, zR
    complex(dp), intent(out) :: Ex(:), Ey(:), Ez(:)

    integer :: i, N
    complex(dp) :: E0
    real(dp) :: khx, khy, khz, px, py, pz

    N = size(obs_x)

!$omp parallel do schedule(static) default(shared) private(i,E0,khx,khy,khz,px,py,pz)
    do i = 1, N
      E0 = gaussian_scalar_field(obs_x(i), obs_y(i), obs_z(i), w0, z0, k0, zR)
      call gaussian_khat_and_polarization( &
          obs_x(i), obs_y(i), obs_z(i), z0, zR, khx, khy, khz, px, py, pz)

      Ex(i) = E0 * px
      Ey(i) = E0 * py
      Ez(i) = E0 * pz
    end do
!$omp end parallel do
  end subroutine calc_gaussian_incident_e_points

  subroutine calc_gaussian_incident_eh_points( &
      obs_x, obs_y, obs_z, w0, z0, k0, zR, eta0, Ex, Ey, Ez, Hx, Hy, Hz)
    real(dp), intent(in) :: obs_x(:), obs_y(:), obs_z(:)
    real(dp), intent(in) :: w0, z0, k0, zR, eta0
    complex(dp), intent(out) :: Ex(:), Ey(:), Ez(:), Hx(:), Hy(:), Hz(:)

    integer :: i, N
    complex(dp) :: E0
    real(dp) :: khx, khy, khz, px, py, pz

    N = size(obs_x)

!$omp parallel do schedule(static) default(shared) private(i,E0,khx,khy,khz,px,py,pz)
    do i = 1, N
      E0 = gaussian_scalar_field(obs_x(i), obs_y(i), obs_z(i), w0, z0, k0, zR)
      call gaussian_khat_and_polarization( &
          obs_x(i), obs_y(i), obs_z(i), z0, zR, khx, khy, khz, px, py, pz)

      Ex(i) = E0 * px
      Ey(i) = E0 * py
      Ez(i) = E0 * pz

      Hx(i) = (khy*Ez(i) - khz*Ey(i)) / eta0
      Hy(i) = (khz*Ex(i) - khx*Ez(i)) / eta0
      Hz(i) = (khx*Ey(i) - khy*Ex(i)) / eta0
    end do
!$omp end parallel do
  end subroutine calc_gaussian_incident_eh_points

  subroutine calc_surface_currents_from_incident_eh( &
      tgt_nx, tgt_ny, tgt_nz, Ex, Ey, Ez, Hx, Hy, Hz, &
      apply_illumination_mask, Jx, Jy, Jz)
    real(dp), intent(in) :: tgt_nx(:), tgt_ny(:), tgt_nz(:)
    complex(dp), intent(in) :: Ex(:), Ey(:), Ez(:), Hx(:), Hy(:), Hz(:)
    logical, intent(in) :: apply_illumination_mask
    complex(dp), intent(out) :: Jx(:), Jy(:), Jz(:)

    integer :: i, N
    real(dp) :: nx, ny, nz, sx, sy, sz
    logical :: illuminated

    N = size(tgt_nx)
    Jx = (0.0_dp, 0.0_dp)
    Jy = (0.0_dp, 0.0_dp)
    Jz = (0.0_dp, 0.0_dp)

!$omp parallel do schedule(static) default(shared) private(i,nx,ny,nz,sx,sy,sz,illuminated)
    do i = 1, N
      nx = tgt_nx(i)
      ny = tgt_ny(i)
      nz = tgt_nz(i)
      illuminated = .true.

      if (apply_illumination_mask) then
        sx = 0.5_dp * real(Ey(i)*conjg(Hz(i)) - Ez(i)*conjg(Hy(i)), dp)
        sy = 0.5_dp * real(Ez(i)*conjg(Hx(i)) - Ex(i)*conjg(Hz(i)), dp)
        sz = 0.5_dp * real(Ex(i)*conjg(Hy(i)) - Ey(i)*conjg(Hx(i)), dp)
        illuminated = (sx*nx + sy*ny + sz*nz) < 0.0_dp
      end if

      if (illuminated) then
        Jx(i) = 2.0_dp * (ny*Hz(i) - nz*Hy(i))
        Jy(i) = 2.0_dp * (nz*Hx(i) - nx*Hz(i))
        Jz(i) = 2.0_dp * (nx*Hy(i) - ny*Hx(i))
      end if
    end do
!$omp end parallel do
  end subroutine calc_surface_currents_from_incident_eh

end module mod_incident
