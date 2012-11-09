! modeule code from http://stackoverflow.com/questions/2813502/reading-input-files-in-fortran
module read_file_m
  use real_m
  implicit none
  contains
  subroutine read_file (UnitNum, FileName, Array)
    integer, intent (in) :: UnitNum
    character (len=*), intent (in) :: FileName
    real(real_t), pointer :: Array(:,:)

    open (unit=UnitNum, file=FileName, status='old', action='read' )
    !backspace (UnitNum)
    block    
      integer :: i, j
      do i=1, size(Array, 2)
        read (UnitNum, *) (Array (i, j), j=1, size(Array, 1))
      end do
    end block
    close (UnitNum)
   end subroutine read_file
end module read_file_m

program test
  use solver_2D_m
  use mpdata_m
  use cyclic_m
  use read_file_m
  implicit none

  if (command_argument_count() /= (9)) then
    print*, "expecting 9 arguments (nx, ny, Cx, Cy, nt, it, f.in, f.out, dec)"
    stop 1
  end if

  block
    type(solver_2D_t) :: slv
    type(mpdata_t), target :: adv
    type(cyclic_t), target :: bcx, bcy
    integer :: nx, ny, nt, it, dec
    real(real_t) :: Cx, Cy

    nx = arg2int(1)
    ny = arg2int(2)
print*, nx, " ", ny
    Cx = arg2real(3)
    Cy = arg2real(4)
    nt = arg2int(5)
    it = arg2int(6)
    dec = arg2int(9)

    call adv%ctor(it)

    block
      class(bcd_t), pointer :: bcx_p, bcy_p
      class(adv_t), pointer :: adv_p
      adv_p => adv
      bcx_p => bcx
      bcy_p => bcy
      call slv%ctor(nx, ny, adv_p, bcx_p, bcy_p)
    end block

    block
      real(real_t), pointer :: tmp(:,:)
      character (len=666) :: fname
      tmp => slv%state() 
      call get_command_argument(7, fname)
      call read_file(666, fname, tmp)

      tmp => slv%Cx() 
      tmp = Cx

      tmp => slv%Cy() 
      tmp = Cy
    end block

    call slv%solve(nt)

    block
      real(real_t), pointer :: tmp(:,:)
      character (len=666) :: fname
      allocate(tmp(nx,ny))
      call get_command_argument(8, fname)
      call read_file(666, fname, tmp)
      if (maxval(abs(slv%state() - tmp)) >= .5 * 10.**(-dec)) stop 1
      deallocate(tmp)
    end block

    call slv%dtor()
    call bcx%dtor()
    call bcy%dtor()
    call adv%dtor()
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
end program
