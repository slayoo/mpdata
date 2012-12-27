program test
  use mpdata_2D_m
  use cyclic_m
  implicit none

  if (command_argument_count() /= (9)) then
    print*, "expecting 9 arguments (nx, ny, Cx, Cy, nt, it, f.in, f.out, dec)"
    stop 1
  end if

  block
    type(mpdata_2D_t) :: slv
    type(cyclic_t), target :: bcx, bcy
    integer :: nx, ny, nt, it, dec
    real(real_t) :: Cx, Cy

    nx = arg2int(1)
    ny = arg2int(2)
    Cx = arg2real(3)
    Cy = arg2real(4)
    nt = arg2int(5)
    it = arg2int(6)
    dec = arg2int(9)

    block
      class(bcd_t), pointer :: bcx_p, bcy_p
      bcx_p => bcx
      bcy_p => bcy
      call slv%ctor(it, bcx_p, bcy_p, nx, ny)
    end block

    block
      real(real_t), pointer :: tmp(:,:)
      character (len=666) :: fname
      tmp => slv%state() 
      call get_command_argument(7, fname)
      call read_file(fname, tmp)

      tmp => slv%courant(0) 
      tmp = Cx

      tmp => slv%courant(1) 
      tmp = Cy
    end block

    call slv%solve(nt)

    block
      logical :: error = .FALSE.
      real(real_t), pointer :: tmp(:,:)
      real :: diff
      character (len=666) :: fname
      allocate(tmp(nx,ny))
      call get_command_argument(8, fname)
      call read_file(fname, tmp)
      diff = maxval(abs(slv%state() - tmp))
      print*, "diff=", diff
      if (diff >= .5 * 10.**(-dec)) error = .TRUE.
      deallocate(tmp)
      if (error) then
        print*, slv%state()
        stop 1
      end if
    end block

    call slv%dtor()
    call bcx%dtor()
    call bcy%dtor()
  end block

  contains

  function arg2int(argno) result(return)
    implicit none
    integer :: argno, return
    character (len=666) :: tmp
    call get_command_argument(argno, tmp)
    read(tmp,*)return
  end function

  function arg2real(argno) result(return)
    implicit none
    integer :: argno
    real(real_t) :: return
    character (len=666) :: tmp
    call get_command_argument(argno, tmp)
    read(tmp,*)return
  end function

  subroutine read_file(fname, array)
    character (len=*), intent (in) :: fname
    real(real_t), pointer :: array(:,:)
    integer :: un
    open(newunit=un, file=fname, status='old', action='read' )
    block    
      integer :: i, j
      do i=1, size(array, 1)
        read(un, *) (array (i, j), j=1, size(array, 2))
      end do
    end block
    close (un)
   end subroutine read_file
end program
