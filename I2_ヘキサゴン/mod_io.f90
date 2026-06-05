module mod_io
  use mod_types
  implicit none
  private
  public :: export_field_map_binary
  public :: export_xy_plane_index
  public :: export_geometry_summary
  public :: export_ipo_summary
  public :: ensure_directory
  public :: create_timestamp_output_dir

contains

  subroutine ensure_directory(dirname)
    character(len=*), intent(in) :: dirname
    integer :: cmdstat, exitstat

    if (len_trim(dirname) <= 0) return
    call execute_command_line( &
        'mkdir -p "' // trim(dirname) // '"', &
        wait=.true., exitstat=exitstat, cmdstat=cmdstat )

    if (cmdstat /= 0 .or. exitstat /= 0) then
      write(*,'(A,A)') 'ERROR: failed to create directory: ', trim(dirname)
      stop 'ensure_directory failed.'
    end if
  end subroutine ensure_directory

  subroutine create_timestamp_output_dir(output_root, output_dir)
    character(len=*), intent(in) :: output_root
    character(len=*), intent(out) :: output_dir

    integer :: values(8)
    character(len=32) :: stamp

    call date_and_time(values=values)
    write(stamp,'(I2.2,I2.2,I2.2,"_",I2.2,I2.2)') &
        mod(values(1), 100), values(2), values(3), values(5), values(6)

    if (len_trim(output_root) > 0) then
      output_dir = trim(output_root) // '/' // trim(stamp)
    else
      output_dir = trim(stamp)
    end if

    call ensure_directory(trim(output_dir))
  end subroutine create_timestamp_output_dir


  subroutine export_field_map_binary(filename, n_pts, obs_x, obs_y, obs_z, e_sq)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n_pts
    real(dp), intent(in) :: obs_x(n_pts), obs_y(n_pts), obs_z(n_pts)
    real(dp), intent(in) :: e_sq(n_pts)
    integer :: u

    open(newunit=u, file=trim(filename), form='unformatted', access='stream', status='replace')
    write(u) obs_x
    write(u) obs_y
    write(u) obs_z
    write(u) e_sq
    close(u)
  end subroutine export_field_map_binary

  subroutine export_xy_plane_index(filename, z_list)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: z_list(:)
    integer :: u, i

    open(newunit=u, file=trim(filename), status='replace', action='write')
    write(u,'(A)') '# index  z[m]'
    do i = 1, size(z_list)
      write(u,'(I4,1X,ES18.10)') i, z_list(i)
    end do
    close(u)
  end subroutine export_xy_plane_index

  subroutine export_geometry_summary( &
      filename, r_cone_in, r_pipe, l_cone, l_pipe, base_z, f_cap, cap_depth, z_cap_vertex, &
      n_total, n_cone, n_pipe, n_cap)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: r_cone_in, r_pipe, l_cone, l_pipe, base_z, f_cap, cap_depth, z_cap_vertex
    integer, intent(in) :: n_total, n_cone, n_pipe, n_cap
    integer :: u

    open(newunit=u, file=trim(filename), status='replace', action='write')
    write(u,'(A)') 'Closed-cap concentrator geometry'
    write(u,'(A,ES18.10)') 'base_z [m]            = ', base_z
    write(u,'(A,ES18.10)') 'r_cone_in [m]         = ', r_cone_in
    write(u,'(A,ES18.10)') 'r_pipe [m]            = ', r_pipe
    write(u,'(A,ES18.10)') 'l_cone [m]            = ', l_cone
    write(u,'(A,ES18.10)') 'l_pipe [m]            = ', l_pipe
    write(u,'(A,ES18.10)') 'f_cap [m]             = ', f_cap
    write(u,'(A,ES18.10)') 'cap_depth [m]         = ', cap_depth
    write(u,'(A,ES18.10)') 'z_cap_vertex [m]      = ', z_cap_vertex
    write(u,'(A,I0)')       'panel_total           = ', n_total
    write(u,'(A,I0)')       'panel_cone            = ', n_cone
    write(u,'(A,I0)')       'panel_pipe            = ', n_pipe
    write(u,'(A,I0)')       'panel_parabolic_cap   = ', n_cap
    close(u)
  end subroutine export_geometry_summary

  subroutine export_ipo_summary(filename, ratios, order_build_times, n_orders_done)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: ratios(:), order_build_times(:)
    integer, intent(in) :: n_orders_done
    integer :: u, i

    open(newunit=u, file=trim(filename), status='replace', action='write')
    write(u,'(A)') '# order  current_ratio  build_time_s'
    do i = 1, n_orders_done
      write(u,'(I4,1X,ES18.10,1X,ES18.10)') i, ratios(i), order_build_times(i)
    end do
    close(u)
  end subroutine export_ipo_summary

end module mod_io
