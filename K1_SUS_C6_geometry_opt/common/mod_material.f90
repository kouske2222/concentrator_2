! exp(-i omega t), inward (metal -> vacuum) normal, SI units.
module mod_material
  use mod_types
  use mod_config, only: sim_config_type
  implicit none
  private
  public :: load_material, surface_current, radiate_panel, zs
  complex(dp) :: zs=(0.0_dp,0.0_dp)
  real(dp) :: eta=376.730313668_dp
contains
  subroutine load_material(filename,cfg)
    character(len=*), intent(in) :: filename
    type(sim_config_type), intent(in) :: cfg
    real(dp) :: conductivity,mu_r
    logical :: pec
    integer :: u,ios
    namelist /material/ conductivity,mu_r,pec
    conductivity=1.37e6_dp; mu_r=1.0_dp; pec=.false.
    open(newunit=u,file=filename,status='old',action='read',iostat=ios)
    if(ios/=0) error stop 'Cannot open material input.'
    read(u,nml=material,iostat=ios); close(u)
    if(ios/=0) error stop 'Missing or invalid /material/ namelist.'
    if(conductivity<=0 .or. mu_r<=0) error stop 'Invalid material constants.'
    eta=cfg%eta0
    zs=(1.0_dp-I_C)*sqrt(PI*cfg%freq*cfg%mu0*mu_r/conductivity)
    if(pec) zs=(0.0_dp,0.0_dp)
    write(*,*) 'Zs [ohm] = ',zs
  end subroutine

  pure function cross(a,b) result(c)
    complex(dp), intent(in) :: a(3),b(3)
    complex(dp) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function

  subroutine surface_current(n,khat,H,J)
    real(dp), intent(in) :: n(3),khat(3)
    complex(dp), intent(in) :: H(3)
    complex(dp), intent(out) :: J(3)
    complex(dp) :: jp(3),s(3),t(3)
    real(dp) :: c,sn
    c=max(0.0_dp,min(1.0_dp,-sum(n*khat)))
    J=(0.0_dp,0.0_dp)
    if(c<=1.0e-12_dp) return
    jp=2.0_dp*cross(cmplx(n,0.0_dp,dp),H)
    s=cross(cmplx(n,0.0_dp,dp),cmplx(khat,0.0_dp,dp))
    sn=sqrt(sum(abs(s)**2))
    if(sn<1.0e-12_dp) then
      J=jp*eta/(eta+zs)
    else
      s=s/sn; t=cross(s,cmplx(n,0.0_dp,dp))
      J=s*sum(s*jp)*eta/(eta+zs*c)+t*sum(t*jp)*(eta*c)/(eta*c+zs)
    end if
  end subroutine

  pure subroutine electric_kernel(k,eta0,rvec,area,J,E,H,propagator)
    real(dp), intent(in) :: k,eta0,rvec(3),area
    complex(dp), intent(in) :: J(3)
    complex(dp), intent(out) :: E(3),H(3)
    complex(dp), intent(in) :: propagator
    real(dp) :: r,u(3),v
    complex(dp) :: phase,a,b,ju
    E=(0.0_dp,0.0_dp); H=E
    r=sqrt(sum(rvec*rvec)); if(r<=tiny(1.0_dp)) return
    v=1.0_dp/r; u=rvec*v; ju=sum(u*J)
    phase=propagator
    a=-I_C*k*v; b=v*v+I_C*v*v*v/k
    E=eta0*phase*(a*(u*ju-J)+b*(3.0_dp*u*ju-J))
    H=phase*(I_C*k*v-v*v)*cross(cmplx(u,0.0_dp,dp),J)
  end subroutine

  subroutine radiate_panel(k,eta0,rvec,area,n,J,E,H,propagator)
    real(dp), intent(in) :: k,eta0,rvec(3),area,n(3)
    complex(dp), intent(in) :: J(3)
    complex(dp), intent(out) :: E(3),H(3)
    complex(dp), optional, intent(in) :: propagator
    complex(dp) :: M(3),Em(3),Hm(3),phase
    if(present(propagator)) then
      phase=propagator
    else
      phase=area*exp(I_C*k*sqrt(sum(rvec*rvec)))/(4.0_dp*PI)
    end if
    call electric_kernel(k,eta0,rvec,area,J,E,H,phase)
    ! Love equivalent currents: J=n x H, M=-n x E=-Zs n x J.
    M=-zs*cross(cmplx(n,0.0_dp,dp),J)
    call electric_kernel(k,eta0,rvec,area,M,Em,Hm,phase)
    E=E-Hm; H=H+Em/(eta0*eta0)
  end subroutine
end module mod_material
