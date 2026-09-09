! Measurement-only entry points; original IPO and radiation routines unchanged.
module timing_bridge
  use iso_c_binding
  use mod_types
  use mod_bridge, only: cfg,mesh
  use mod_incident, only: build_initial_po_current
  use mod_geometry, only: radius_at_z
  use mod_operator, only: field_from_current
  implicit none
contains
  subroutine bm_initial(j,areas) bind(C)
    complex(c_double_complex), intent(out) :: j(3,mesh%Q,6)
    real(c_double), intent(out) :: areas(mesh%Q,6)
    real(dp) :: norms(6)
    integer :: p
    call build_initial_po_current(mesh,cfg,j,norms)
    do p=1,6
      areas(:,p)=mesh%area(p,:)
    end do
  end subroutine
  subroutine bm_fields(j,e,h) bind(C)
    complex(c_double_complex), intent(in) :: j(3,mesh%Q,6)
    complex(c_double_complex), intent(out) :: e(3,189),h(3,189)
    real(dp) :: x(189),y(189),z(189),zz,rr,phi
    integer :: i,iz,jj
    i=0
    do iz=0,20
      zz=cfg%base_z+(0.01_dp+0.98_dp*real(iz,dp)/20.0_dp)*(cfg%l_cone+cfg%l_pipe)
      rr=0.7_dp*radius_at_z(zz,cfg)*cos(PI/6.0_dp)
      do jj=0,8
        i=i+1;z(i)=zz;x(i)=0;y(i)=0
        if(jj>0) then
          phi=2*PI*real(jj-1,dp)/8
          x(i)=rr*cos(phi);y(i)=rr*sin(phi)
        end if
      end do
    end do
    call field_from_current(mesh,cfg,j,x,y,z,e,h,.false.)
  end subroutine
end module
