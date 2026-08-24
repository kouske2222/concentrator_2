module mod_config
  use mod_types
  implicit none
  private
  public :: sim_config_type, load_config, load_validation_config, update_derived_config
  public :: estimate_problem, write_preflight

  type :: sim_config_type
    real(dp) :: c0, mu0, eps0, eta0
    real(dp) :: freq, lam0, k0
    real(dp) :: e0, w0_src, z0_src, zR_src
    real(dp) :: beam_offset_x, beam_offset_y, incidence_theta_deg, incidence_phi_deg, polarization_deg
    logical :: use_plane_wave
    integer :: M
    real(dp) :: d_in, d_out, l_cone, l_pipe, base_z
    logical :: diameter_across_vertices
    integer :: n_face, n_z_cone, n_z_pipe
    integer :: max_reflection_order
    integer :: run_mode                  ! 1=full only, 2=C6 modal only, 3=postprocess compare
    real(dp) :: stop_ratio
    integer :: visibility_samples
    real(dp) :: wall_clearance_fraction
    integer :: nx_xz, nz_xz, nxy, naxis
    logical :: use_dense_full, run_full, compute_fields, save_order_exit_fields
    logical :: write_mode_fields
    real(dp) :: max_dense_gib, estimated_pair_rate
    integer :: mode_policy                 ! 0=all, 1=m=1,5, 2=automatic
    real(dp) :: mode_tolerance
    character(len=256) :: output_dir
  end type sim_config_type

contains

  subroutine set_defaults(cfg)
    type(sim_config_type), intent(out) :: cfg
    cfg%c0=2.99792458e8_dp; cfg%mu0=4.0e-7_dp*PI
    cfg%eps0=1.0_dp/(cfg%mu0*cfg%c0*cfg%c0); cfg%eta0=sqrt(cfg%mu0/cfg%eps0)
    cfg%freq=3.0e9_dp; cfg%e0=1.0_dp; cfg%w0_src=0.080_dp; cfg%z0_src=-0.25_dp
    cfg%beam_offset_x=0.0_dp; cfg%beam_offset_y=0.0_dp
    cfg%incidence_theta_deg=0.0_dp; cfg%incidence_phi_deg=0.0_dp; cfg%polarization_deg=0.0_dp
    cfg%use_plane_wave=.false.; cfg%M=6
    cfg%d_in=0.180_dp; cfg%d_out=0.065_dp; cfg%l_cone=0.160_dp; cfg%l_pipe=0.100_dp
    cfg%base_z=0.0_dp; cfg%diameter_across_vertices=.true.
    cfg%n_face=3; cfg%n_z_cone=6; cfg%n_z_pipe=4
    cfg%max_reflection_order=4; cfg%run_mode=3; cfg%stop_ratio=1.0e-4_dp; cfg%visibility_samples=7
    cfg%wall_clearance_fraction=0.04_dp
    cfg%nx_xz=51; cfg%nz_xz=81; cfg%nxy=41; cfg%naxis=161
    cfg%use_dense_full=.true.; cfg%run_full=.true.; cfg%compute_fields=.true.
    cfg%save_order_exit_fields=.true.; cfg%write_mode_fields=.true.; cfg%max_dense_gib=2.0_dp
    cfg%estimated_pair_rate=2.0e7_dp; cfg%mode_policy=0; cfg%mode_tolerance=1.0e-8_dp
    cfg%output_dir='results/validation'
    call update_derived_config(cfg)
  end subroutine set_defaults

  subroutine load_validation_config(cfg)
    type(sim_config_type), intent(out) :: cfg
    call set_defaults(cfg)
  end subroutine load_validation_config

  subroutine load_config(filename,cfg)
    character(len=*), intent(in) :: filename
    type(sim_config_type), intent(out) :: cfg
    integer :: unit,ios
    real(dp) :: freq,e0,w0_src,z0_src,d_in,d_out,l_cone,l_pipe,base_z
    real(dp) :: beam_offset_x,beam_offset_y,incidence_theta_deg,incidence_phi_deg,polarization_deg
    real(dp) :: stop_ratio,wall_clearance_fraction,max_dense_gib,estimated_pair_rate,mode_tolerance
    integer :: M,n_face,n_z_cone,n_z_pipe,max_reflection_order,run_mode,visibility_samples
    integer :: nx_xz,nz_xz,nxy,naxis,mode_policy
    logical :: use_plane_wave,diameter_across_vertices,use_dense_full,run_full,compute_fields
    logical :: save_order_exit_fields,write_mode_fields
    character(len=256) :: output_dir
    namelist /simulation/ freq,e0,w0_src,z0_src,beam_offset_x,beam_offset_y,incidence_theta_deg, &
      incidence_phi_deg,polarization_deg,use_plane_wave,M,d_in,d_out,l_cone,l_pipe,base_z, &
      diameter_across_vertices,n_face,n_z_cone,n_z_pipe,max_reflection_order,run_mode,stop_ratio, &
      visibility_samples,wall_clearance_fraction,nx_xz,nz_xz,nxy,naxis,use_dense_full,run_full, &
      compute_fields,save_order_exit_fields,write_mode_fields,max_dense_gib,estimated_pair_rate,mode_policy, &
      mode_tolerance,output_dir

    call set_defaults(cfg)
    freq=cfg%freq; e0=cfg%e0; w0_src=cfg%w0_src; z0_src=cfg%z0_src
    beam_offset_x=cfg%beam_offset_x; beam_offset_y=cfg%beam_offset_y
    incidence_theta_deg=cfg%incidence_theta_deg; incidence_phi_deg=cfg%incidence_phi_deg
    polarization_deg=cfg%polarization_deg
    use_plane_wave=cfg%use_plane_wave; M=cfg%M; d_in=cfg%d_in; d_out=cfg%d_out
    l_cone=cfg%l_cone; l_pipe=cfg%l_pipe; base_z=cfg%base_z
    diameter_across_vertices=cfg%diameter_across_vertices; n_face=cfg%n_face
    n_z_cone=cfg%n_z_cone; n_z_pipe=cfg%n_z_pipe; max_reflection_order=cfg%max_reflection_order
    run_mode=cfg%run_mode
    stop_ratio=cfg%stop_ratio; visibility_samples=cfg%visibility_samples
    wall_clearance_fraction=cfg%wall_clearance_fraction; nx_xz=cfg%nx_xz; nz_xz=cfg%nz_xz
    nxy=cfg%nxy; naxis=cfg%naxis; use_dense_full=cfg%use_dense_full; run_full=cfg%run_full
    compute_fields=cfg%compute_fields; save_order_exit_fields=cfg%save_order_exit_fields
    write_mode_fields=cfg%write_mode_fields
    max_dense_gib=cfg%max_dense_gib; estimated_pair_rate=cfg%estimated_pair_rate
    mode_policy=cfg%mode_policy; mode_tolerance=cfg%mode_tolerance; output_dir=cfg%output_dir
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios/=0) then
      write(*,'(A)') 'Cannot open configuration namelist: '//trim(filename)
      stop 1
    end if
    read(unit,nml=simulation,iostat=ios); close(unit)
    if (ios/=0) then
      write(*,'(A)') 'Invalid /simulation/ namelist: '//trim(filename)
      stop 1
    end if
    cfg%freq=freq; cfg%e0=e0; cfg%w0_src=w0_src; cfg%z0_src=z0_src
    cfg%beam_offset_x=beam_offset_x; cfg%beam_offset_y=beam_offset_y
    cfg%incidence_theta_deg=incidence_theta_deg; cfg%incidence_phi_deg=incidence_phi_deg
    cfg%polarization_deg=polarization_deg
    cfg%use_plane_wave=use_plane_wave; cfg%M=M; cfg%d_in=d_in; cfg%d_out=d_out
    cfg%l_cone=l_cone; cfg%l_pipe=l_pipe; cfg%base_z=base_z
    cfg%diameter_across_vertices=diameter_across_vertices; cfg%n_face=n_face
    cfg%n_z_cone=n_z_cone; cfg%n_z_pipe=n_z_pipe; cfg%max_reflection_order=max_reflection_order
    cfg%run_mode=run_mode
    cfg%stop_ratio=stop_ratio; cfg%visibility_samples=visibility_samples
    cfg%wall_clearance_fraction=wall_clearance_fraction; cfg%nx_xz=nx_xz; cfg%nz_xz=nz_xz
    cfg%nxy=nxy; cfg%naxis=naxis; cfg%use_dense_full=use_dense_full; cfg%run_full=run_full
    cfg%compute_fields=compute_fields; cfg%save_order_exit_fields=save_order_exit_fields
    cfg%write_mode_fields=write_mode_fields
    cfg%max_dense_gib=max_dense_gib; cfg%estimated_pair_rate=estimated_pair_rate
    cfg%mode_policy=mode_policy; cfg%mode_tolerance=mode_tolerance; cfg%output_dir=output_dir
    call update_derived_config(cfg)
    if (cfg%M/=6) then
      write(*,'(A,I0,A)') 'A regular hexagonal concentrator requires M=6 (received M=',cfg%M,').'
      stop 1
    end if
  end subroutine load_config

  subroutine update_derived_config(cfg)
    type(sim_config_type), intent(inout) :: cfg
    cfg%lam0=cfg%c0/cfg%freq; cfg%k0=2.0_dp*PI/cfg%lam0
    cfg%zR_src=PI*cfg%w0_src**2/cfg%lam0
  end subroutine update_derived_config

  subroutine estimate_problem(cfg,Q,N,dense_gib,c6_gib,pairs_full,pairs_c6,seconds_full)
    type(sim_config_type), intent(in) :: cfg
    integer, intent(out) :: Q,N
    real(dp), intent(out) :: dense_gib,c6_gib,seconds_full
    integer(int64), intent(out) :: pairs_full,pairs_c6
    Q=2*cfg%n_face*(cfg%n_z_cone+cfg%n_z_pipe); N=cfg%M*Q
    pairs_full=int(N,int64)*int(N,int64); pairs_c6=int(cfg%M,int64)*int(Q,int64)*int(Q,int64)
    dense_gib=16.0_dp*real((3_int64*N)*(3_int64*N),dp)/1024.0_dp**3
    c6_gib=16.0_dp*real((3_int64*Q)*(3_int64*Q)*cfg%M,dp)/1024.0_dp**3
    seconds_full=real(pairs_full,dp)/max(cfg%estimated_pair_rate,1.0_dp)
  end subroutine estimate_problem

  subroutine write_preflight(cfg,filename)
    type(sim_config_type), intent(in) :: cfg
    character(len=*), intent(in) :: filename
    integer :: Q,N,u
    integer(int64) :: pf,pc
    real(dp) :: gf,gc,seconds
    call estimate_problem(cfg,Q,N,gf,gc,pf,pc,seconds)
    open(newunit=u,file=filename,status='replace',action='write')
    write(u,'(A,ES24.16)') 'frequency_hz=',cfg%freq
    write(u,'(A,ES24.16)') 'wavelength_m=',cfg%lam0
    write(u,'(A,I0)') 'panels_per_sector=',Q; write(u,'(A,I0)') 'total_panels=',N
    write(u,'(A,I0)') 'full_pair_interactions_per_order=',pf
    write(u,'(A,I0)') 'c6_pair_interactions_per_order=',pc
    write(u,'(A,ES24.16)') 'full_dense_gib=',gf; write(u,'(A,ES24.16)') 'c6_blocks_gib=',gc
    write(u,'(A,I0)') 'max_reflection_order=',cfg%max_reflection_order
    write(u,'(A,ES24.16)') 'estimated_full_build_seconds=',seconds
    write(u,'(A,ES24.16)') 'estimated_full_ipo_seconds=',seconds*real(cfg%max_reflection_order,dp)
    write(u,'(A,L1)') 'full_dense_safe=',gf<=cfg%max_dense_gib
    close(u)
  end subroutine write_preflight
end module mod_config
