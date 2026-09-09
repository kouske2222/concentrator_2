program test_surface
  use mod_types
  use mod_config
  use mod_material
  implicit none
  type(sim_config_type) :: cfg
  real(dp) :: theta,c,s,n(3),k(3)
  complex(dp) :: h(3),j(3),expected(3),rte,rtm
  integer :: i
  call load_config('config_smoke.nml',cfg)
  call load_material('config_smoke.nml',cfg)
  if(real(zs)<=0 .or. aimag(zs)>=0) error stop 'Passive impedance/phasor sign failed'
  n=[0.0_dp,0.0_dp,1.0_dp]
  do i=0,89
    theta=real(i,dp)*PI/180.0_dp; c=cos(theta); s=sin(theta)
    k=[s,0.0_dp,-c]
    h=cmplx([c,0.0_dp,s]/cfg%eta0,0.0_dp,dp)
    call surface_current(n,k,h,j)
    expected=(0.0_dp,0.0_dp); expected(2)=2.0_dp*c/(cfg%eta0+zs*c)
    if(maxval(abs(j-expected))>1.0e-12_dp) error stop 'TE current mismatch'
    h=(0.0_dp,0.0_dp); h(2)=-1.0_dp/cfg%eta0
    call surface_current(n,k,h,j)
    expected=(0.0_dp,0.0_dp); expected(1)=2.0_dp*c/(cfg%eta0*c+zs)
    if(maxval(abs(j-expected))>1.0e-12_dp) error stop 'TM current mismatch'
    rte=(zs*c-cfg%eta0)/(zs*c+cfg%eta0)
    rtm=(zs-cfg%eta0*c)/(zs+cfg%eta0*c)
    if(abs(rte)>1.0_dp .or. abs(rtm)>1.0_dp) error stop 'Nonpassive Fresnel reflection'
  end do
  print *, 'TE/TM plane-wave currents and passive reflection: PASS'
end program
