program run_94ghz_sweep
  use mod_types
  use mod_config
  implicit none
  type(sim_config_type) :: cfg
  character(len=256) :: config_path,csv_path
  real(dp), parameter :: frequencies(5)=[90.0e9_dp,92.0e9_dp,94.0e9_dp,96.0e9_dp,98.0e9_dp]
  real(dp), parameter :: polarizations(3)=[0.0_dp,45.0_dp,90.0_dp]
  real(dp), parameter :: offsets_mm(3)=[0.0_dp,1.0_dp,2.0_dp]
  real(dp), parameter :: angles_deg(3)=[0.0_dp,0.5_dp,1.0_dp]
  real(dp) :: full_gib,c8_gib,seconds,input_freq
  integer(int64) :: full_pairs,c8_pairs
  integer :: i,j,k,l,Q,N,u

  config_path='config_94ghz.nml'
  if (command_argument_count()>=1) call get_command_argument(1,config_path)
  call load_config(trim(config_path),cfg)
  input_freq=cfg%freq
  call execute_command_line('mkdir -p '//trim(cfg%output_dir))
  call write_preflight(cfg,trim(cfg%output_dir)//'/preflight_94ghz.txt')
  csv_path=trim(cfg%output_dir)//'/sweep_preflight.csv'
  open(newunit=u,file=trim(csv_path),status='replace',action='write')
  write(u,'(A)') 'frequency_hz,polarization_deg,offset_mm,incidence_deg,panels_per_sector,total_panels,'// &
                 'full_pairs_per_order,c8_pairs_per_order,full_dense_gib,c8_blocks_gib,'// &
                 'operator_reusable_for_condition,full_dense_safe'
  do i=1,size(frequencies)
    cfg%freq=frequencies(i); call update_derived_config(cfg)
    call estimate_problem(cfg,Q,N,full_gib,c8_gib,full_pairs,c8_pairs,seconds)
    do j=1,size(polarizations); do k=1,size(offsets_mm); do l=1,size(angles_deg)
      write(u,'(ES16.8,3(",",F8.3),2(",",I0),2(",",I0),2(",",ES16.8),",",L1,",",L1)') &
        cfg%freq,polarizations(j),offsets_mm(k),angles_deg(l),Q,N,full_pairs,c8_pairs, &
        full_gib,c8_gib,.true.,full_gib<=cfg%max_dense_gib
    end do; end do; end do
  end do
  close(u)
  cfg%freq=input_freq; call update_derived_config(cfg)
  write(*,'(A)') '94 GHz sweep preflight completed: '//trim(csv_path)
  write(*,'(A,I0)') 'Conditions listed = ',size(frequencies)*size(polarizations)*size(offsets_mm)*size(angles_deg)
  write(*,'(A,ES12.4)') 'C8 panel pairs per reflection/order = ',real(c8_pairs,dp)
  if (c8_gib>cfg%max_dense_gib) then
    write(*,'(A)') 'Dense C8 storage is unsafe. Use matrix-free GPU/FMM or an admissible-block H-matrix/ACA.'
    write(*,'(A)') 'The operator is reusable across polarization, offset, and angle at fixed frequency.'
    write(*,'(A)') 'A frequency change requires rebuilding the Green-function propagation operator.'
  end if
end program run_94ghz_sweep
