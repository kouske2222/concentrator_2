module mod_checkpoint
  use mod_types
  implicit none
  private
  public :: ensure_checkpoint_dir
  public :: save_ipo_checkpoint_meta, load_ipo_checkpoint_meta
  public :: save_ipo_checkpoint_order, load_ipo_checkpoint_order

contains

  subroutine ensure_checkpoint_dir(checkpoint_dir)
    character(len=*), intent(in) :: checkpoint_dir
    integer :: cmdstat, exitstat

    if (len_trim(checkpoint_dir) <= 0) return
    call execute_command_line( &
        'mkdir -p "' // trim(checkpoint_dir) // '"', &
        wait=.true., exitstat=exitstat, cmdstat=cmdstat )
  end subroutine ensure_checkpoint_dir

  function ipo_checkpoint_meta_filename(checkpoint_dir) result(fname)
    character(len=*), intent(in) :: checkpoint_dir
    character(len=512) :: fname
    write(fname,'(A,"/ipo_meta.bin")') trim(checkpoint_dir)
  end function ipo_checkpoint_meta_filename

  function ipo_checkpoint_order_filename(checkpoint_dir, order) result(fname)
    character(len=*), intent(in) :: checkpoint_dir
    integer, intent(in) :: order
    character(len=512) :: fname
    write(fname,'(A,"/ipo_order_",I6.6,".bin")') trim(checkpoint_dir), order
  end function ipo_checkpoint_order_filename

  subroutine save_ipo_checkpoint_meta( &
      checkpoint_dir, N, max_order, last_order, ratios, order_build_times)
    character(len=*), intent(in) :: checkpoint_dir
    integer, intent(in) :: N, max_order, last_order
    real(dp), intent(in) :: ratios(:), order_build_times(:)

    integer :: u
    character(len=8) :: magic
    character(len=512) :: fname

    magic = 'IPOMETA1'
    fname = ipo_checkpoint_meta_filename(checkpoint_dir)

    open(newunit=u, file=trim(fname), form='unformatted', access='stream', status='replace')
    write(u) magic
    write(u) N
    write(u) max_order
    write(u) last_order
    write(u) ratios(1:max_order)
    write(u) order_build_times(1:max_order)
    close(u)
  end subroutine save_ipo_checkpoint_meta

  subroutine load_ipo_checkpoint_meta( &
      checkpoint_dir, N_expected, max_order_expected, last_order, ratios, order_build_times, ok)
    character(len=*), intent(in) :: checkpoint_dir
    integer, intent(in) :: N_expected, max_order_expected
    integer, intent(out) :: last_order
    real(dp), intent(out) :: ratios(:), order_build_times(:)
    logical, intent(out) :: ok

    integer :: u, N_file, max_order_file
    character(len=8) :: magic
    character(len=512) :: fname

    fname = ipo_checkpoint_meta_filename(checkpoint_dir)
    inquire(file=trim(fname), exist=ok)

    if (.not. ok) then
      last_order = 0
      return
    end if

    open(newunit=u, file=trim(fname), form='unformatted', access='stream', status='old')
    read(u) magic
    read(u) N_file
    read(u) max_order_file
    read(u) last_order

    if (magic /= 'IPOMETA1') stop 'Invalid IPO checkpoint meta magic.'
    if (N_file /= N_expected) stop 'IPO checkpoint mesh size mismatch. Delete the checkpoint directory.'
    if (max_order_file /= max_order_expected) stop 'IPO checkpoint max_order mismatch. Delete the checkpoint directory.'

    read(u) ratios(1:max_order_expected)
    read(u) order_build_times(1:max_order_expected)
    close(u)

    if (last_order < 0 .or. last_order > max_order_expected) then
      stop 'Invalid IPO checkpoint last_order.'
    end if
  end subroutine load_ipo_checkpoint_meta

  subroutine save_ipo_checkpoint_order(checkpoint_dir, order, Jx, Jy, Jz)
    character(len=*), intent(in) :: checkpoint_dir
    integer, intent(in) :: order
    complex(dp), intent(in) :: Jx(:), Jy(:), Jz(:)

    integer :: u, N
    character(len=8) :: magic
    character(len=512) :: fname

    magic = 'IPOORD01'
    N = size(Jx)
    fname = ipo_checkpoint_order_filename(checkpoint_dir, order)

    open(newunit=u, file=trim(fname), form='unformatted', access='stream', status='replace')
    write(u) magic
    write(u) order
    write(u) N
    write(u) Jx
    write(u) Jy
    write(u) Jz
    close(u)
  end subroutine save_ipo_checkpoint_order

  subroutine load_ipo_checkpoint_order(checkpoint_dir, order, Jx, Jy, Jz)
    character(len=*), intent(in) :: checkpoint_dir
    integer, intent(in) :: order
    complex(dp), intent(out) :: Jx(:), Jy(:), Jz(:)

    integer :: u, order_file, N_file
    character(len=8) :: magic
    character(len=512) :: fname

    fname = ipo_checkpoint_order_filename(checkpoint_dir, order)

    open(newunit=u, file=trim(fname), form='unformatted', access='stream', status='old')
    read(u) magic
    read(u) order_file
    read(u) N_file

    if (magic /= 'IPOORD01') stop 'Invalid IPO checkpoint order magic.'
    if (order_file /= order) stop 'IPO checkpoint order mismatch.'
    if (N_file /= size(Jx)) stop 'IPO checkpoint order size mismatch.'

    read(u) Jx
    read(u) Jy
    read(u) Jz
    close(u)
  end subroutine load_ipo_checkpoint_order

end module mod_checkpoint
