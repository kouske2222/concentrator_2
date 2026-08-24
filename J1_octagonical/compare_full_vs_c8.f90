program compare_full_vs_c8
  use omp_lib, only: omp_get_max_threads
  use mod_types
  use mod_config
  use mod_geometry
  use mod_incident
  use mod_operator
  use mod_full3d
  use mod_c8_modal
  implicit none

  type(sim_config_type) :: cfg
  type(panel_mesh_type) :: mesh
  complex(dp), allocatable :: J0(:,:,:)
  complex(dp), allocatable :: full_orders(:,:,:,:),modal_orders(:,:,:,:)
  complex(dp), allocatable :: Jfull(:,:,:),Jmodal(:,:,:)
  real(dp), allocatable :: rfull(:),rmodal(:),tfull(:),tmodal(:)
  real(dp), allocatable :: order_errors(:)
  real(dp) :: initial_modes(8),modal_initial_modes(8),modal_final_modes(8),fast_final_modes(8)
  real(dp) :: full_build_s,block_build_s,circ_error,current_error,fast_error,total_t0,total_t1
  real(dp) :: xz_error,axis_error,exit_error,exit_power_full,exit_power_modal
  real(dp) :: dense_gib,c8_gib,estimated_seconds
  real(dp) :: full_total_s,modal_total_s,full_step_sum_s,modal_step_sum_s
  integer(int64) :: pairs_full,pairs_c8,full_elements
  integer :: nf,nm,Q_est,N_est,full_done,modal_done
  character(len=256) :: config_path
  character(len=512) :: full_binary_path,modal_binary_path

  call cpu_time(total_t0)
  config_path='config_validation.nml'
  if (command_argument_count()>=1) call get_command_argument(1,config_path)
  call load_config(trim(config_path),cfg)
  call execute_command_line('mkdir -p '//trim(cfg%output_dir))
  call write_preflight(cfg,trim(cfg%output_dir)//'/preflight.txt')
  call estimate_problem(cfg,Q_est,N_est,dense_gib,c8_gib,pairs_full,pairs_c8,estimated_seconds)
  write(*,'(A,F12.3)') 'Estimated full dense storage [GiB] = ',dense_gib
  write(*,'(A)') 'Estimated C8 block storage [GiB] =        0.000 (matrix-free)'
  write(*,'(A,ES12.4)') 'Full panel pairs per order         = ',real(pairs_full,dp)
  write(*,'(A,I0)') 'Run mode                          = ',cfg%run_mode
  write(*,'(A,I0)') 'OpenMP max threads                = ',omp_get_max_threads()
  call build_c8_mesh(cfg,mesh)
  full_binary_path=trim(cfg%output_dir)//'/J_full_output.bin'
  modal_binary_path=trim(cfg%output_dir)//'/J_modal_output.bin'

  full_build_s=0.0_dp; block_build_s=0.0_dp; circ_error=0.0_dp
  current_error=0.0_dp; fast_error=0.0_dp
  xz_error=0.0_dp; axis_error=0.0_dp; exit_error=0.0_dp
  exit_power_full=0.0_dp; exit_power_modal=0.0_dp
  full_total_s=0.0_dp; modal_total_s=0.0_dp; full_step_sum_s=0.0_dp; modal_step_sum_s=0.0_dp
  full_elements=0_int64; nf=1; nm=1
  full_done=0; modal_done=0
  initial_modes=0.0_dp; modal_initial_modes=0.0_dp; modal_final_modes=0.0_dp; fast_final_modes=0.0_dp

  select case (cfg%run_mode)
  case (1)
    write(*,'(A)') 'Run mode 1: Full 3D matrix-free IPO only.'
    allocate(J0(3,mesh%Q,mesh%M))
    call build_initial_po_current(mesh,cfg,J0,initial_modes)
    call run_full_ipo_matrix_free(mesh,cfg,J0,cfg%max_reflection_order,cfg%stop_ratio, &
                                  full_orders,rfull,tfull,nf)
    allocate(Jfull(3,mesh%Q,mesh%M))
    Jfull=sum(full_orders(:,:,:,1:nf),dim=4)
    call write_current_binary(trim(full_binary_path),Jfull,mesh%Q,mesh%M)
    call cpu_time(total_t1)
    call write_runtime_summary(cfg,'full',nf,rfull,tfull,total_t1-total_t0)
    call write_metrics(cfg,mesh,circ_error,current_error,fast_error,xz_error,axis_error,exit_error, &
                       exit_power_full,exit_power_modal,full_build_s,block_build_s,total_t1-total_t0, &
                       full_elements,0_int64,rfull,nf)
    write(*,'(A)') 'Full current binary written: '//trim(full_binary_path)

  case (2)
    write(*,'(A)') 'Run mode 2: C8 modal matrix-free IPO only.'
    allocate(J0(3,mesh%Q,mesh%M))
    call build_initial_po_current(mesh,cfg,J0,initial_modes)
    call run_c8_modal_ipo(mesh,cfg,J0,cfg%max_reflection_order,cfg%stop_ratio, &
                          modal_orders,rmodal,tmodal,modal_initial_modes,modal_final_modes,nm)
    allocate(Jmodal(3,mesh%Q,mesh%M))
    Jmodal=sum(modal_orders(:,:,:,1:nm),dim=4)
    fast_final_modes=modal_final_modes
    call write_current_binary(trim(modal_binary_path),Jmodal,mesh%Q,mesh%M)
    call write_mode_csv(cfg,initial_modes,modal_final_modes,fast_final_modes)
    call cpu_time(total_t1)
    call write_runtime_summary(cfg,'modal',nm,rmodal,tmodal,total_t1-total_t0)
    call write_metrics(cfg,mesh,circ_error,current_error,fast_error,xz_error,axis_error,exit_error, &
                       exit_power_full,exit_power_modal,full_build_s,block_build_s,total_t1-total_t0, &
                       full_elements,0_int64,rmodal,nm)
    write(*,'(A)') 'Modal current binary written: '//trim(modal_binary_path)

  case (3)
    write(*,'(A)') 'Run mode 3: postprocess comparison from saved current binaries.'
    allocate(Jfull(3,mesh%Q,mesh%M),Jmodal(3,mesh%Q,mesh%M))
    call cpu_time(full_build_s)
    call read_current_binary(trim(full_binary_path),Jfull,mesh%Q,mesh%M)
    call cpu_time(block_build_s)
    full_build_s=block_build_s-full_build_s
    call read_current_binary(trim(modal_binary_path),Jmodal,mesh%Q,mesh%M)
    call cpu_time(total_t1)
    block_build_s=total_t1-block_build_s
    current_error=array_relative_error(Jmodal,Jfull)
    fast_error=current_error
    allocate(rfull(1),rmodal(1),tfull(1),tmodal(1),order_errors(1))
    rfull(1)=0.0_dp; rmodal(1)=0.0_dp; tfull(1)=full_build_s; tmodal(1)=block_build_s
    order_errors(1)=current_error
    if (cfg%compute_fields) then
      call compare_and_write_fields(mesh,cfg,Jfull,Jmodal,xz_error,axis_error,exit_error, &
                                    exit_power_full,exit_power_modal)
    end if
    call write_postprocess_csv(cfg,current_error,full_build_s,block_build_s)
    call read_runtime_binary(trim(cfg%output_dir)//'/runtime_full.bin',full_done,full_total_s,full_step_sum_s)
    call read_runtime_binary(trim(cfg%output_dir)//'/runtime_modal.bin',modal_done,modal_total_s,modal_step_sum_s)
    call write_runtime_comparison_csv(cfg,full_done,modal_done,full_total_s,modal_total_s, &
                                      full_step_sum_s,modal_step_sum_s)
    call write_current_csv(cfg,mesh,Jfull,Jmodal)
    call cpu_time(total_t1)
    call write_metrics(cfg,mesh,circ_error,current_error,fast_error,xz_error,axis_error,exit_error, &
                       exit_power_full,exit_power_modal,full_build_s,block_build_s,total_t1-total_t0, &
                       full_elements,0_int64,rfull,1)

  case default
    write(*,'(A,I0)') 'Invalid run_mode. Use 1=full, 2=C8 modal, or 3=postprocess. run_mode = ',cfg%run_mode
    stop 1
  end select

  write(*,'(A)') 'C8 PO/IPO comparison completed.'
  write(*,'(A,I0,A,I0)') 'Panels: total=',mesh%N,', per sector=',mesh%Q
  if (cfg%run_mode==3) then
    write(*,'(A,ES12.4)') 'Cumulative full-modal error   = ',current_error
    write(*,'(A,ES12.4)') 'XZ field relative error       = ',xz_error
  end if
  write(*,'(A)') 'Results: '//trim(cfg%output_dir)

contains

  real(dp) function array_relative_error(a,b) result(value)
    complex(dp), intent(in) :: a(:,:,:),b(:,:,:)
    value=sqrt(sum(abs(a-b)**2))/(sqrt(sum(abs(b)**2))+1.0e-300_dp)
  end function array_relative_error

  subroutine write_current_binary(filename,J,Q,M)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: Q,M
    complex(dp), intent(in) :: J(3,Q,M)
    character(len=8) :: magic
    integer :: u,ios
    magic='JCURR01 '
    open(newunit=u,file=trim(filename),status='replace',form='unformatted',action='write',iostat=ios)
    if (ios/=0) then
      write(*,'(A)') 'Cannot open current binary for writing: '//trim(filename)
      stop 1
    end if
    write(u,iostat=ios) magic
    if (ios/=0) stop 'Failed to write current binary magic.'
    write(u,iostat=ios) Q
    if (ios/=0) stop 'Failed to write current binary Q.'
    write(u,iostat=ios) M
    if (ios/=0) stop 'Failed to write current binary M.'
    write(u,iostat=ios) J
    if (ios/=0) stop 'Failed to write current binary data.'
    close(u,iostat=ios)
    if (ios/=0) stop 'Failed to close current binary after writing.'
  end subroutine write_current_binary

  subroutine read_current_binary(filename,J,Q,M)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: Q,M
    complex(dp), intent(out) :: J(3,Q,M)
    character(len=8) :: magic
    integer :: u,ios,Q_file,M_file
    logical :: exists
    inquire(file=trim(filename),exist=exists)
    if (.not.exists) then
      write(*,'(A)') 'Current binary does not exist: '//trim(filename)
      stop 1
    end if
    open(newunit=u,file=trim(filename),status='old',form='unformatted',action='read',iostat=ios)
    if (ios/=0) then
      write(*,'(A)') 'Cannot open current binary for reading: '//trim(filename)
      stop 1
    end if
    read(u,iostat=ios) magic
    if (ios/=0) stop 'Failed to read current binary magic.'
    if (magic/='JCURR01 ') stop 'Invalid current binary magic.'
    read(u,iostat=ios) Q_file
    if (ios/=0) stop 'Failed to read current binary Q.'
    read(u,iostat=ios) M_file
    if (ios/=0) stop 'Failed to read current binary M.'
    if (Q_file/=Q) then
      write(*,'(A,I0,A,I0)') 'Current binary Q mismatch: file=',Q_file,', expected=',Q
      stop 1
    end if
    if (M_file/=M) then
      write(*,'(A,I0,A,I0)') 'Current binary M mismatch: file=',M_file,', expected=',M
      stop 1
    end if
    read(u,iostat=ios) J
    if (ios/=0) stop 'Failed to read current binary data.'
    close(u,iostat=ios)
    if (ios/=0) stop 'Failed to close current binary after reading.'
  end subroutine read_current_binary

  subroutine write_postprocess_csv(conf,current_error,full_read_s,modal_read_s)
    type(sim_config_type), intent(in) :: conf
    real(dp), intent(in) :: current_error,full_read_s,modal_read_s
    integer :: u
    open(newunit=u,file=trim(conf%output_dir)//'/postprocess_summary.csv',status='replace',action='write')
    write(u,'(A)') 'current_rel_error,full_current_read_s,modal_current_read_s'
    write(u,'(3(ES24.16,:,","))') current_error,full_read_s,modal_read_s
    close(u)
  end subroutine write_postprocess_csv

  subroutine write_runtime_summary(conf,label,n_done,ratios,step_times,total_s)
    type(sim_config_type), intent(in) :: conf
    character(len=*), intent(in) :: label
    integer, intent(in) :: n_done
    real(dp), intent(in) :: ratios(:),step_times(:),total_s
    character(len=512) :: csv_name,bin_name
    real(dp) :: step_sum_s,avg_step_s,last_ratio
    integer :: u,i

    step_sum_s=0.0_dp
    if (n_done>=2) step_sum_s=sum(step_times(2:n_done))
    if (n_done>=2) then
      avg_step_s=step_sum_s/real(n_done-1,dp)
    else
      avg_step_s=0.0_dp
    end if
    last_ratio=0.0_dp
    if (n_done>=1) last_ratio=ratios(n_done)

    csv_name=trim(conf%output_dir)//'/runtime_'//trim(label)//'.csv'
    open(newunit=u,file=trim(csv_name),status='replace',action='write')
    write(u,'(A)') 'label,run_mode,max_reflection_order,n_orders_done,total_cpu_s,ipo_step_sum_cpu_s,avg_ipo_step_cpu_s,last_ratio'
    write(u,'(A,",",I0,",",I0,",",I0,4(",",ES24.16))') trim(label),conf%run_mode, &
      conf%max_reflection_order,n_done,total_s,step_sum_s,avg_step_s,last_ratio
    write(u,'(A)') 'order,ratio,step_cpu_s'
    do i=1,n_done
      write(u,'(I0,2(",",ES24.16))') i-1,ratios(i),step_times(i)
    end do
    close(u)

    bin_name=trim(conf%output_dir)//'/runtime_'//trim(label)//'.bin'
    call write_runtime_binary(trim(bin_name),n_done,total_s,step_sum_s)
  end subroutine write_runtime_summary

  subroutine write_runtime_binary(filename,n_done,total_s,step_sum_s)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n_done
    real(dp), intent(in) :: total_s,step_sum_s
    character(len=8) :: magic
    integer :: u,ios
    magic='RTIME01 '
    open(newunit=u,file=trim(filename),status='replace',form='unformatted',action='write',iostat=ios)
    if (ios/=0) then
      write(*,'(A)') 'Cannot open runtime binary for writing: '//trim(filename)
      stop 1
    end if
    write(u,iostat=ios) magic
    if (ios/=0) stop 'Failed to write runtime binary magic.'
    write(u,iostat=ios) n_done
    if (ios/=0) stop 'Failed to write runtime binary n_done.'
    write(u,iostat=ios) total_s
    if (ios/=0) stop 'Failed to write runtime binary total_s.'
    write(u,iostat=ios) step_sum_s
    if (ios/=0) stop 'Failed to write runtime binary step_sum_s.'
    close(u,iostat=ios)
    if (ios/=0) stop 'Failed to close runtime binary after writing.'
  end subroutine write_runtime_binary

  subroutine read_runtime_binary(filename,n_done,total_s,step_sum_s)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n_done
    real(dp), intent(out) :: total_s,step_sum_s
    character(len=8) :: magic
    integer :: u,ios
    logical :: exists
    inquire(file=trim(filename),exist=exists)
    if (.not.exists) then
      write(*,'(A)') 'Runtime binary does not exist: '//trim(filename)
      stop 1
    end if
    open(newunit=u,file=trim(filename),status='old',form='unformatted',action='read',iostat=ios)
    if (ios/=0) then
      write(*,'(A)') 'Cannot open runtime binary for reading: '//trim(filename)
      stop 1
    end if
    read(u,iostat=ios) magic
    if (ios/=0) stop 'Failed to read runtime binary magic.'
    if (magic/='RTIME01 ') stop 'Invalid runtime binary magic.'
    read(u,iostat=ios) n_done
    if (ios/=0) stop 'Failed to read runtime binary n_done.'
    read(u,iostat=ios) total_s
    if (ios/=0) stop 'Failed to read runtime binary total_s.'
    read(u,iostat=ios) step_sum_s
    if (ios/=0) stop 'Failed to read runtime binary step_sum_s.'
    close(u,iostat=ios)
    if (ios/=0) stop 'Failed to close runtime binary after reading.'
  end subroutine read_runtime_binary

  subroutine write_runtime_comparison_csv(conf,full_done,modal_done,full_total_s,modal_total_s, &
                                          full_step_sum_s,modal_step_sum_s)
    type(sim_config_type), intent(in) :: conf
    integer, intent(in) :: full_done,modal_done
    real(dp), intent(in) :: full_total_s,modal_total_s,full_step_sum_s,modal_step_sum_s
    real(dp) :: total_speedup,ipo_speedup
    integer :: u
    total_speedup=full_total_s/(modal_total_s+1.0e-300_dp)
    ipo_speedup=full_step_sum_s/(modal_step_sum_s+1.0e-300_dp)
    open(newunit=u,file=trim(conf%output_dir)//'/runtime_comparison.csv',status='replace',action='write')
    write(u,'(A)') 'quantity,full,modal,full_over_modal'
    write(u,'(A,",",I0,",",I0,",",ES24.16)') 'n_orders_done',full_done,modal_done, &
      real(full_done,dp)/(real(modal_done,dp)+1.0e-300_dp)
    write(u,'(A,3(",",ES24.16))') 'total_cpu_s',full_total_s,modal_total_s,total_speedup
    write(u,'(A,3(",",ES24.16))') 'ipo_step_sum_cpu_s',full_step_sum_s,modal_step_sum_s,ipo_speedup
    close(u)
  end subroutine write_runtime_comparison_csv

  subroutine compare_and_write_fields(msh,conf,Jf,Jm,err_xz,err_axis,err_exit,pf,pm)
    type(panel_mesh_type), intent(in) :: msh
    type(sim_config_type), intent(in) :: conf
    complex(dp), intent(in) :: Jf(3,msh%Q,msh%M),Jm(3,msh%Q,msh%M)
    real(dp), intent(out) :: err_xz,err_axis,err_exit,pf,pm
    real(dp), allocatable :: x(:),y(:),z(:)
    complex(dp), allocatable :: Ef(:,:),Hf(:,:),Em(:,:),Hm(:,:)
    logical, allocatable :: mask(:)
    integer :: i,ix,iz,n,plane
    real(dp) :: rin,dx,dz,zplanes(3),zexit,cell_area
    character(len=512) :: filef,filem
    rin=0.5_dp*conf%d_in

    n=conf%nx_xz*conf%nz_xz
    allocate(x(n),y(n),z(n),mask(n),Ef(3,n),Hf(3,n),Em(3,n),Hm(3,n))
    dx=2.0_dp*rin/real(conf%nx_xz-1,dp); dz=(conf%l_cone+conf%l_pipe)/real(conf%nz_xz-1,dp)
    i=0
    do iz=1,conf%nz_xz; do ix=1,conf%nx_xz
      i=i+1; x(i)=-rin+real(ix-1,dp)*dx; y(i)=0.0_dp; z(i)=conf%base_z+real(iz-1,dp)*dz
      mask(i)=inside_octagon(x(i),y(i),z(i),conf,conf%wall_clearance_fraction)
    end do; end do
    call field_from_current(msh,conf,Jf,x,y,z,Ef,Hf,.true.)
    call field_from_current(msh,conf,Jm,x,y,z,Em,Hm,.true.)
    call zero_masked(Ef,Hf,mask); call zero_masked(Em,Hm,mask)
    err_xz=field_error(Em,Ef,mask)
    call write_field_csv(trim(conf%output_dir)//'/xz_full.csv',x,y,z,Ef,mask)
    call write_field_csv(trim(conf%output_dir)//'/xz_modal.csv',x,y,z,Em,mask)
    deallocate(x,y,z,mask,Ef,Hf,Em,Hm)

    n=conf%naxis; allocate(x(n),y(n),z(n),mask(n),Ef(3,n),Hf(3,n),Em(3,n),Hm(3,n))
    do i=1,n
      x(i)=0.0_dp; y(i)=0.0_dp; z(i)=conf%base_z+(conf%l_cone+conf%l_pipe)*real(i-1,dp)/real(n-1,dp); mask(i)=.true.
    end do
    call field_from_current(msh,conf,Jf,x,y,z,Ef,Hf,.true.)
    call field_from_current(msh,conf,Jm,x,y,z,Em,Hm,.true.)
    err_axis=field_error(Em,Ef,mask)
    call write_field_csv(trim(conf%output_dir)//'/axis_full.csv',x,y,z,Ef,mask)
    call write_field_csv(trim(conf%output_dir)//'/axis_modal.csv',x,y,z,Em,mask)
    deallocate(x,y,z,mask,Ef,Hf,Em,Hm)

    zplanes=[conf%base_z+0.5_dp*conf%l_cone,conf%base_z+conf%l_cone+0.25_dp*conf%l_pipe, &
             conf%base_z+conf%l_cone+0.75_dp*conf%l_pipe]
    do plane=1,3
      call make_xy_grid(conf,zplanes(plane),x,y,z,mask,dx)
      n=size(x); allocate(Ef(3,n),Hf(3,n),Em(3,n),Hm(3,n))
      call field_from_current(msh,conf,Jf,x,y,z,Ef,Hf,.true.)
      call field_from_current(msh,conf,Jm,x,y,z,Em,Hm,.true.)
      call zero_masked(Ef,Hf,mask); call zero_masked(Em,Hm,mask)
      write(filef,'(A,"/xy",I0,"_full.csv")') trim(conf%output_dir),plane
      write(filem,'(A,"/xy",I0,"_modal.csv")') trim(conf%output_dir),plane
      call write_field_csv(trim(filef),x,y,z,Ef,mask); call write_field_csv(trim(filem),x,y,z,Em,mask)
      deallocate(x,y,z,mask,Ef,Hf,Em,Hm)
    end do

    zexit=conf%base_z+conf%l_cone+conf%l_pipe-0.25_dp*min(conf%l_cone/real(conf%n_z_cone,dp), &
                                                        conf%l_pipe/real(conf%n_z_pipe,dp))
    call make_xy_grid(conf,zexit,x,y,z,mask,dx); n=size(x)
    allocate(Ef(3,n),Hf(3,n),Em(3,n),Hm(3,n))
    call field_from_current(msh,conf,Jf,x,y,z,Ef,Hf,.true.)
    call field_from_current(msh,conf,Jm,x,y,z,Em,Hm,.true.)
    call zero_masked(Ef,Hf,mask); call zero_masked(Em,Hm,mask)
    err_exit=field_error(Em,Ef,mask); cell_area=dx*dx
    pf=exit_power(Ef,Hf,mask,cell_area); pm=exit_power(Em,Hm,mask,cell_area)
    call write_field_csv(trim(conf%output_dir)//'/exit_full.csv',x,y,z,Ef,mask)
    call write_field_csv(trim(conf%output_dir)//'/exit_modal.csv',x,y,z,Em,mask)
    deallocate(x,y,z,mask,Ef,Hf,Em,Hm)
  end subroutine compare_and_write_fields

  subroutine make_xy_grid(conf,zplane,x,y,z,mask,spacing)
    type(sim_config_type), intent(in) :: conf
    real(dp), intent(in) :: zplane
    real(dp), allocatable, intent(out) :: x(:),y(:),z(:)
    logical, allocatable, intent(out) :: mask(:)
    real(dp), intent(out) :: spacing
    integer :: n,i,ix,iy
    real(dp) :: radius
    radius=radius_at_z(zplane,conf); n=conf%nxy*conf%nxy
    allocate(x(n),y(n),z(n),mask(n)); spacing=2.0_dp*radius/real(conf%nxy-1,dp); i=0
    do iy=1,conf%nxy; do ix=1,conf%nxy
      i=i+1; x(i)=-radius+real(ix-1,dp)*spacing; y(i)=-radius+real(iy-1,dp)*spacing; z(i)=zplane
      mask(i)=inside_octagon(x(i),y(i),z(i),conf,conf%wall_clearance_fraction)
    end do; end do
  end subroutine make_xy_grid

  subroutine zero_masked(E,H,mask)
    complex(dp), intent(inout) :: E(:,:),H(:,:)
    logical, intent(in) :: mask(:)
    integer :: i
    do i=1,size(mask)
      if (.not.mask(i)) then; E(:,i)=(0.0_dp,0.0_dp); H(:,i)=(0.0_dp,0.0_dp); end if
    end do
  end subroutine zero_masked

  real(dp) function field_error(a,b,mask) result(value)
    complex(dp), intent(in) :: a(:,:),b(:,:)
    logical, intent(in) :: mask(:)
    integer :: i
    real(dp) :: num,den
    num=0.0_dp; den=0.0_dp
    do i=1,size(mask)
      if (mask(i)) then; num=num+sum(abs(a(:,i)-b(:,i))**2); den=den+sum(abs(b(:,i))**2); end if
    end do
    value=sqrt(num)/(sqrt(den)+1.0e-300_dp)
  end function field_error

  real(dp) function exit_power(E,H,mask,area) result(power)
    complex(dp), intent(in) :: E(:,:),H(:,:)
    logical, intent(in) :: mask(:)
    real(dp), intent(in) :: area
    integer :: i
    power=0.0_dp
    do i=1,size(mask)
      if (mask(i)) power=power+0.5_dp*real(E(1,i)*conjg(H(2,i))-E(2,i)*conjg(H(1,i)),dp)*area
    end do
  end function exit_power

  subroutine write_field_csv(filename,x,y,z,E,mask)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: x(:),y(:),z(:)
    complex(dp), intent(in) :: E(:,:)
    logical, intent(in) :: mask(:)
    integer :: u,i
    open(newunit=u,file=filename,status='replace',action='write')
    write(u,'(A)') 'x_m,y_m,z_m,mask,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im,E_mag,Ex_phase_rad'
    do i=1,size(x)
      write(u,'(3(ES24.16,","),I1,",",8(ES24.16,:,","))') x(i),y(i),z(i),merge(1,0,mask(i)), &
        real(E(1,i)),aimag(E(1,i)),real(E(2,i)),aimag(E(2,i)),real(E(3,i)),aimag(E(3,i)), &
        sqrt(sum(abs(E(:,i))**2)),atan2(aimag(E(1,i)),real(E(1,i)))
    end do
    close(u)
  end subroutine write_field_csv

  subroutine write_mode_csv(conf,initial,final,fastfinal)
    type(sim_config_type), intent(in) :: conf
    real(dp), intent(in) :: initial(:),final(:),fastfinal(:)
    integer :: u,m
    open(newunit=u,file=trim(conf%output_dir)//'/mode_norms.csv',status='replace')
    write(u,'(A)') 'mode,initial_norm,final_norm,fast_final_norm'
    do m=1,size(initial); write(u,'(I0,3(",",ES24.16))') m-1,initial(m),final(m),fastfinal(m); end do
    close(u)
  end subroutine write_mode_csv

  subroutine write_current_csv(conf,msh,Jf,Jm)
    type(sim_config_type), intent(in) :: conf
    type(panel_mesh_type), intent(in) :: msh
    complex(dp), intent(in) :: Jf(3,msh%Q,msh%M),Jm(3,msh%Q,msh%M)
    integer :: u,p,q
    open(newunit=u,file=trim(conf%output_dir)//'/surface_current.csv',status='replace')
    write(u,'(A)') 'sector,panel,x_m,y_m,z_m,Jfull_mag,Jmodal_mag'
    do p=1,msh%M; do q=1,msh%Q
      write(u,'(I0,",",I0,5(",",ES24.16))') p-1,q,msh%x(p,q),msh%y(p,q),msh%z(p,q), &
        sqrt(sum(abs(Jf(:,q,p))**2)),sqrt(sum(abs(Jm(:,q,p))**2))
    end do; end do; close(u)
  end subroutine write_current_csv

  subroutine write_metrics(conf,msh,ecirc,ej,efast,exz,eaxis,eexit,pf,pm,tbuild,tblock,ttotal, &
                           full_count,c8_count,ratios,nf_l)
    type(sim_config_type), intent(in) :: conf
    type(panel_mesh_type), intent(in) :: msh
    real(dp), intent(in) :: ecirc,ej,efast,exz,eaxis,eexit,pf,pm,tbuild,tblock,ttotal,ratios(:)
    integer(int64), intent(in) :: full_count,c8_count
    integer, intent(in) :: nf_l
    integer :: u
    real(dp) :: full_mib,block_mib
    full_mib=16.0_dp*real(full_count,dp)/1024.0_dp**2
    block_mib=16.0_dp*real(c8_count,dp)/1024.0_dp**2
    open(newunit=u,file=trim(conf%output_dir)//'/metrics.txt',status='replace')
    write(u,'(A,ES24.16)') 'frequency_hz=',conf%freq
    write(u,'(A,I0)') 'run_mode=',conf%run_mode
    write(u,'(A,I0)') 'M=',conf%M; write(u,'(A,I0)') 'panels_per_sector=',msh%Q; write(u,'(A,I0)') 'total_panels=',msh%N
    write(u,'(A,ES24.16)') 'block_circulant_relative_error=',ecirc
    write(u,'(A,ES24.16)') 'cumulative_full_modal_relative_error=',ej
    write(u,'(A,ES24.16)') 'cumulative_full_fast_relative_error=',efast
    write(u,'(A,ES24.16)') 'xz_field_relative_error=',exz
    write(u,'(A,ES24.16)') 'axis_field_relative_error=',eaxis
    write(u,'(A,ES24.16)') 'exit_field_relative_error=',eexit
    write(u,'(A,ES24.16)') 'exit_power_full_w=',pf; write(u,'(A,ES24.16)') 'exit_power_modal_w=',pm
    write(u,'(A,ES24.16)') 'full_operator_build_s=',tbuild; write(u,'(A,ES24.16)') 'c8_blocks_build_s=',tblock
    if (conf%run_mode==3) then
      write(u,'(A,ES24.16)') 'full_current_read_s=',tbuild
      write(u,'(A,ES24.16)') 'modal_current_read_s=',tblock
    end if
    write(u,'(A,ES24.16)') 'total_elapsed_s=',ttotal
    write(u,'(A,I0)') 'full_stored_complex_elements=',full_count
    write(u,'(A,I0)') 'c8_stored_complex_elements=',c8_count
    write(u,'(A,ES24.16)') 'full_operator_mib=',full_mib; write(u,'(A,ES24.16)') 'c8_blocks_mib=',block_mib
    write(u,'(A,ES24.16)') 'last_current_ratio=',ratios(nf_l)
    write(u,'(A,L1)') 'ipo_converged=',ratios(nf_l)<conf%stop_ratio
    write(u,'(A,F8.3)') 'ideal_speedup_all_modes=',real(conf%M,dp)
    write(u,'(A,F8.3)') 'ideal_speedup_modes_1_7=',real(conf%M*conf%M,dp)/2.0_dp
    close(u)
  end subroutine write_metrics

end program compare_full_vs_c8
