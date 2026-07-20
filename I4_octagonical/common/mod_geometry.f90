module mod_geometry
  use mod_types
  use mod_config, only: sim_config_type
  implicit none
  private
  public :: panel_mesh_type, build_c8_mesh, radius_at_z, inside_octagon
  public :: rotation_z, global_to_local, local_to_global

  type :: panel_mesh_type
    integer :: M, Q, N
    real(dp), allocatable :: x(:,:), y(:,:), z(:,:)
    real(dp), allocatable :: nx(:,:), ny(:,:), nz(:,:)
    real(dp), allocatable :: area(:,:), panel_len(:,:)
    real(dp), allocatable :: vertices(:,:,:,:)
  end type panel_mesh_type

contains

  pure real(dp) function circumradius(diameter, across_vertices) result(radius)
    real(dp), intent(in) :: diameter
    logical, intent(in) :: across_vertices
    if (across_vertices) then
      radius = 0.5_dp * diameter
    else
      radius = 0.5_dp * diameter / cos(PI/8.0_dp)
    end if
  end function circumradius

  pure real(dp) function radius_at_z(z, cfg) result(radius)
    real(dp), intent(in) :: z
    type(sim_config_type), intent(in) :: cfg
    real(dp) :: rin, rout, t
    rin = circumradius(cfg%d_in, cfg%diameter_across_vertices)
    rout = circumradius(cfg%d_out, cfg%diameter_across_vertices)
    if (z <= cfg%base_z + cfg%l_cone) then
      t = min(max((z-cfg%base_z)/cfg%l_cone, 0.0_dp), 1.0_dp)
      radius = rin + t*(rout-rin)
    else
      radius = rout
    end if
  end function radius_at_z

  pure subroutine rotation_z(phi, R)
    real(dp), intent(in) :: phi
    real(dp), intent(out) :: R(3,3)
    real(dp) :: c, s
    c = cos(phi); s = sin(phi)
    R = 0.0_dp
    R(1,1)=c; R(1,2)=-s
    R(2,1)=s; R(2,2)= c
    R(3,3)=1.0_dp
  end subroutine rotation_z

  pure subroutine local_to_global(phi, vl, vg)
    real(dp), intent(in) :: phi
    complex(dp), intent(in) :: vl(3)
    complex(dp), intent(out) :: vg(3)
    real(dp) :: c, s
    c=cos(phi); s=sin(phi)
    vg(1)=c*vl(1)-s*vl(2)
    vg(2)=s*vl(1)+c*vl(2)
    vg(3)=vl(3)
  end subroutine local_to_global

  pure subroutine global_to_local(phi, vg, vl)
    real(dp), intent(in) :: phi
    complex(dp), intent(in) :: vg(3)
    complex(dp), intent(out) :: vl(3)
    real(dp) :: c, s
    c=cos(phi); s=sin(phi)
    vl(1)= c*vg(1)+s*vg(2)
    vl(2)=-s*vg(1)+c*vg(2)
    vl(3)=vg(3)
  end subroutine global_to_local

  pure subroutine rotate_real(phi, vin, vout)
    real(dp), intent(in) :: phi, vin(3)
    real(dp), intent(out) :: vout(3)
    real(dp) :: c, s
    c=cos(phi); s=sin(phi)
    vout(1)=c*vin(1)-s*vin(2)
    vout(2)=s*vin(1)+c*vin(2)
    vout(3)=vin(3)
  end subroutine rotate_real

  subroutine build_c8_mesh(cfg, mesh)
    type(sim_config_type), intent(in) :: cfg
    type(panel_mesh_type), intent(out) :: mesh
    integer :: iz, iu, tri, q, p, nztot
    real(dp), allocatable :: zlev(:), grid(:,:,:)
    real(dp) :: u, r, phi, v0(3), v1(3), v2(3)
    real(dp) :: ea(3), eb(3), cr(3), normc, base_vertices(3,3)
    real(dp), allocatable :: bx(:),by(:),bz(:),bnx(:),bny(:),bnz(:),ba(:),bplen(:)
    real(dp), allocatable :: bv(:,:,:)

    nztot = cfg%n_z_cone + cfg%n_z_pipe
    mesh%M = cfg%M
    mesh%Q = 2*cfg%n_face*nztot
    mesh%N = mesh%M*mesh%Q
    allocate(zlev(nztot+1), grid(nztot+1,cfg%n_face+1,3))
    do iz=0,cfg%n_z_cone
      zlev(iz+1)=cfg%base_z + cfg%l_cone*real(iz,dp)/real(cfg%n_z_cone,dp)
    end do
    do iz=1,cfg%n_z_pipe
      zlev(cfg%n_z_cone+iz+1)=cfg%base_z+cfg%l_cone+cfg%l_pipe*real(iz,dp)/real(cfg%n_z_pipe,dp)
    end do
    do iz=1,nztot+1
      r=radius_at_z(zlev(iz),cfg)
      do iu=0,cfg%n_face
        u=real(iu,dp)/real(cfg%n_face,dp)
        grid(iz,iu+1,1)=r*((1.0_dp-u)*cos(-PI/8.0_dp)+u*cos(PI/8.0_dp))
        grid(iz,iu+1,2)=r*((1.0_dp-u)*sin(-PI/8.0_dp)+u*sin(PI/8.0_dp))
        grid(iz,iu+1,3)=zlev(iz)
      end do
    end do
    allocate(bx(mesh%Q),by(mesh%Q),bz(mesh%Q),bnx(mesh%Q),bny(mesh%Q),bnz(mesh%Q))
    allocate(ba(mesh%Q),bplen(mesh%Q),bv(mesh%Q,3,3))
    q=0
    do iz=1,nztot
      do iu=1,cfg%n_face
        do tri=1,2
          q=q+1
          if (tri==1) then
            base_vertices(1,:)=grid(iz,iu,:)
            base_vertices(2,:)=grid(iz,iu+1,:)
            base_vertices(3,:)=grid(iz+1,iu+1,:)
          else
            base_vertices(1,:)=grid(iz,iu,:)
            base_vertices(2,:)=grid(iz+1,iu+1,:)
            base_vertices(3,:)=grid(iz+1,iu,:)
          end if
          bv(q,:,:)=base_vertices
          bx(q)=sum(base_vertices(:,1))/3.0_dp
          by(q)=sum(base_vertices(:,2))/3.0_dp
          bz(q)=sum(base_vertices(:,3))/3.0_dp
          ea=base_vertices(2,:)-base_vertices(1,:)
          eb=base_vertices(3,:)-base_vertices(1,:)
          cr(1)=ea(2)*eb(3)-ea(3)*eb(2)
          cr(2)=ea(3)*eb(1)-ea(1)*eb(3)
          cr(3)=ea(1)*eb(2)-ea(2)*eb(1)
          normc=sqrt(sum(cr*cr))
          ba(q)=0.5_dp*normc
          bplen(q)=sqrt(ba(q)/PI)
          bnx(q)=-cr(1)/normc; bny(q)=-cr(2)/normc; bnz(q)=-cr(3)/normc
        end do
      end do
    end do
    allocate(mesh%x(mesh%M,mesh%Q),mesh%y(mesh%M,mesh%Q),mesh%z(mesh%M,mesh%Q))
    allocate(mesh%nx(mesh%M,mesh%Q),mesh%ny(mesh%M,mesh%Q),mesh%nz(mesh%M,mesh%Q))
    allocate(mesh%area(mesh%M,mesh%Q),mesh%panel_len(mesh%M,mesh%Q))
    allocate(mesh%vertices(mesh%M,mesh%Q,3,3))
    do p=0,mesh%M-1
      phi=2.0_dp*PI*real(p,dp)/real(mesh%M,dp)
      do q=1,mesh%Q
        v0=[bx(q),by(q),bz(q)]; call rotate_real(phi,v0,v1)
        mesh%x(p+1,q)=v1(1); mesh%y(p+1,q)=v1(2); mesh%z(p+1,q)=v1(3)
        v0=[bnx(q),bny(q),bnz(q)]; call rotate_real(phi,v0,v1)
        mesh%nx(p+1,q)=v1(1); mesh%ny(p+1,q)=v1(2); mesh%nz(p+1,q)=v1(3)
        mesh%area(p+1,q)=ba(q); mesh%panel_len(p+1,q)=bplen(q)
        do tri=1,3
          v0=bv(q,tri,:)
          call rotate_real(phi,v0,v2)
          mesh%vertices(p+1,q,tri,:)=v2
        end do
      end do
    end do
    deallocate(zlev,grid,bx,by,bz,bnx,bny,bnz,ba,bplen,bv)
  end subroutine build_c8_mesh

  pure logical function inside_octagon(x,y,z,cfg,clearance_fraction) result(ok)
    real(dp), intent(in) :: x,y,z,clearance_fraction
    type(sim_config_type), intent(in) :: cfg
    integer :: p
    real(dp) :: phi, apothem, projection
    ok=.false.
    if (z<cfg%base_z .or. z>cfg%base_z+cfg%l_cone+cfg%l_pipe) return
    apothem=radius_at_z(z,cfg)*cos(PI/8.0_dp)*(1.0_dp-clearance_fraction)
    do p=0,cfg%M-1
      phi=2.0_dp*PI*real(p,dp)/real(cfg%M,dp)
      projection=x*cos(phi)+y*sin(phi)
      if (projection>apothem+2.0e-12_dp) return
    end do
    ok=.true.
  end function inside_octagon

end module mod_geometry
