program worker
  use mod_types
  use mod_config
  use mod_geometry
  use mod_incident
  use mod_material
  use mod_operator_old
  use mod_modal_old
  implicit none
  type(sim_config_type) :: cfg
  type(panel_mesh_type) :: mesh
  character(len=2048) :: filename,action
  complex(dp), allocatable :: J(:,:,:),Jn(:,:,:),Jh(:,:,:),Jhn(:,:,:)
  integer, parameter :: np=189
  complex(dp) :: E(3,np),H(3,np),Ei(3,np),Hi(3),temp(3)
  real(dp) :: x(np),y(np),z(np),khat(3),norms(6),rr,phi,zz
  integer :: u,i,jj,p,iz,modes(6),ios
  call get_command_argument(1,filename); call get_command_argument(2,action)
  call load_config(trim(filename),cfg); call load_material(trim(filename),cfg)
  call build_c6_mesh(cfg,mesh)
  allocate(J(3,mesh%Q,6),Jn(3,mesh%Q,6),Jh(3,mesh%Q,6),Jhn(3,mesh%Q,6))
  modes=[0,1,2,3,4,5]
  if(trim(action)=='init') then
    call build_initial_po_current(mesh,cfg,J,norms)
  else
    open(newunit=u,file='input.bin',status='old',access='stream',form='unformatted',iostat=ios)
    if(ios/=0) error stop 'Cannot read input.bin'
    read(u,iostat=ios) J; close(u)
    if(ios/=0) error stop 'Truncated input.bin'
    select case(trim(action))
    case('step')
      call decompose_current_all_modes(J,Jh,norms)
      call apply_c6_modal_operator_matrix_free(mesh,cfg,modes,Jh,Jhn)
      call reconstruct_active_modes(modes,Jhn,Jn)
      J=Jn
    case('full')
      call apply_full_operator_matrix_free(mesh,cfg,J,Jn)
      J=Jn
    case('fields')
      ! Postprocess supplied cumulative current without advancing.
    case default
      error stop 'Action must be init, step, full, or fields.'
    end select
  end if
  ! 21 axial levels, each with axis + 8 points at 70% local apothem.
  i=0
  do iz=0,20
    zz=cfg%base_z+(0.01_dp+0.98_dp*real(iz,dp)/20.0_dp)*(cfg%l_cone+cfg%l_pipe)
    rr=0.7_dp*radius_at_z(zz,cfg)*cos(PI/6.0_dp)
    do jj=0,8
      i=i+1; z(i)=zz; x(i)=0; y(i)=0
      if(jj>0) then
        phi=2.0_dp*PI*real(jj-1,dp)/8.0_dp
        x(i)=rr*cos(phi); y(i)=rr*sin(phi)
      end if
      call incident_eh_point(x(i),y(i),z(i),cfg,Ei(:,i),Hi,khat)
    end do
  end do
  call field_from_current(mesh,cfg,J,x,y,z,E,H,.false.)
  open(newunit=u,file='current.bin',status='replace',access='stream',form='unformatted')
  write(u) J; close(u)
  open(newunit=u,file='electric.bin',status='replace',access='stream',form='unformatted')
  write(u) E; close(u)
  open(newunit=u,file='incident.bin',status='replace',access='stream',form='unformatted')
  write(u) Ei; close(u)
  open(newunit=u,file='areas.bin',status='replace',access='stream',form='unformatted')
  do p=1,6
    write(u) mesh%area(p,:)
  end do
  close(u)
  open(newunit=u,file='probes.bin',status='replace',access='stream',form='unformatted')
  do i=1,np
    write(u) x(i),y(i),z(i)
  end do
  close(u)
end program
