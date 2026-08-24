module mod_incident
  use mod_types
  use mod_config, only: sim_config_type
  use mod_geometry, only: panel_mesh_type, global_to_local
  implicit none
  private
  public :: incident_eh_point, build_initial_po_current

contains

  pure subroutine cross_real_complex(a,b,c)
    real(dp), intent(in) :: a(3)
    complex(dp), intent(in) :: b(3)
    complex(dp), intent(out) :: c(3)
    c(1)=a(2)*b(3)-a(3)*b(2); c(2)=a(3)*b(1)-a(1)*b(3); c(3)=a(1)*b(2)-a(2)*b(1)
  end subroutine cross_real_complex

  pure subroutine incident_frame(cfg,axis,e1)
    type(sim_config_type), intent(in) :: cfg
    real(dp), intent(out) :: axis(3),e1(3)
    real(dp) :: eperp(3)
    real(dp) :: theta,phi,pol,c,s
    theta=cfg%incidence_theta_deg*PI/180.0_dp; phi=cfg%incidence_phi_deg*PI/180.0_dp
    axis=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
    e1=[cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)]
    eperp=[-sin(phi),cos(phi),0.0_dp]
    pol=cfg%polarization_deg*PI/180.0_dp; c=cos(pol); s=sin(pol)
    e1=c*e1+s*eperp
  end subroutine incident_frame

  pure subroutine incident_eh_point(x,y,z,cfg,E,H,khat)
    real(dp), intent(in) :: x,y,z
    type(sim_config_type), intent(in) :: cfg
    complex(dp), intent(out) :: E(3),H(3)
    real(dp), intent(out) :: khat(3)
    real(dp) :: axis(3),pol(3),origin(3),rel(3),trans(3),dz,r2,wz,psi,invR,Rz
    real(dp) :: amp,phase,kn,dot0,pn
    complex(dp) :: scalar
    call incident_frame(cfg,axis,pol)
    origin=[cfg%beam_offset_x,cfg%beam_offset_y,cfg%z0_src]
    rel=[x,y,z]-origin; dz=sum(rel*axis); trans=rel-dz*axis; r2=sum(trans*trans)
    if (cfg%use_plane_wave) then
      khat=axis; scalar=cfg%e0*exp(I_C*cfg%k0*dz)
    else
      wz=cfg%w0_src*sqrt(1.0_dp+(dz/cfg%zR_src)**2); psi=atan(dz/cfg%zR_src)
      if (abs(dz)<1.0e-14_dp) then
        invR=0.0_dp
      else
        Rz=dz*(1.0_dp+(cfg%zR_src/dz)**2); invR=1.0_dp/Rz
      end if
      amp=cfg%e0*(cfg%w0_src/wz)*exp(-r2/(wz*wz))
      phase=cfg%k0*dz+0.5_dp*cfg%k0*r2*invR-psi; scalar=amp*exp(I_C*phase)
      khat=axis+trans*invR; kn=sqrt(sum(khat*khat)); khat=khat/kn
    end if
    dot0=sum(pol*khat); pol=pol-dot0*khat; pn=sqrt(sum(pol*pol)); pol=pol/pn
    E=scalar*cmplx(pol,0.0_dp,kind=dp); call cross_real_complex(khat,E,H); H=H/cfg%eta0
  end subroutine incident_eh_point

  subroutine build_initial_po_current(mesh,cfg,Jlocal,mode_norms)
    type(panel_mesh_type), intent(in) :: mesh
    type(sim_config_type), intent(in) :: cfg
    complex(dp), intent(out) :: Jlocal(3,mesh%Q,mesh%M)
    real(dp), intent(out) :: mode_norms(mesh%M)
    integer :: p,q,m
    real(dp) :: phi,n(3),khat(3),sn
    complex(dp) :: E(3),H(3),Jg(3),Jl(3),phase_m
    complex(dp), allocatable :: Jhat(:,:,:)
    allocate(Jhat(3,mesh%Q,mesh%M)); Jlocal=(0.0_dp,0.0_dp)
    do p=0,mesh%M-1
      phi=2.0_dp*PI*real(p,dp)/real(mesh%M,dp)
      do q=1,mesh%Q
        call incident_eh_point(mesh%x(p+1,q),mesh%y(p+1,q),mesh%z(p+1,q),cfg,E,H,khat)
        n=[mesh%nx(p+1,q),mesh%ny(p+1,q),mesh%nz(p+1,q)]; sn=sum(khat*n)
        if (sn<0.0_dp) then
          Jg(1)=2.0_dp*(n(2)*H(3)-n(3)*H(2))
          Jg(2)=2.0_dp*(n(3)*H(1)-n(1)*H(3))
          Jg(3)=2.0_dp*(n(1)*H(2)-n(2)*H(1))
          call global_to_local(phi,Jg,Jl); Jlocal(:,q,p+1)=Jl
        end if
      end do
    end do
    Jhat=(0.0_dp,0.0_dp)
    do m=0,mesh%M-1
      do p=0,mesh%M-1
        phase_m=exp(-I_C*2.0_dp*PI*real(m*p,dp)/real(mesh%M,dp))
        Jhat(:,:,m+1)=Jhat(:,:,m+1)+Jlocal(:,:,p+1)*phase_m
      end do
      mode_norms(m+1)=sqrt(sum(abs(Jhat(:,:,m+1))**2))
    end do
    deallocate(Jhat)
  end subroutine build_initial_po_current
end module mod_incident
