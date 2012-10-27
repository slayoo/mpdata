program test
  use solver_2D_m
  use mpdata_m
  use cyclic_m
  implicit none

  type(solver_2D_t) :: slv
  type(mpdata_t), target :: adv
  type(cyclic_t), target :: bcx, bcy
  integer :: nx = 3, ny = 3, nt = 1

  call adv%ctor()

  block
    class(bcd_t), pointer :: bcx_p, bcy_p
    class(adv_t), pointer :: adv_p
    adv_p => adv
    bcx_p => bcx
    bcy_p => bcy
    call slv%ctor(nx, ny, adv_p, bcx_p, bcy_p)
  end block

  call init()
  print*, ""
  print '(3F6.3)', slv%state()
  call slv%solve(1)
  print*, ""
  print '(3F6.3)', slv%state()
  call slv%solve(1)
  print*, ""
  print '(3F6.3)', slv%state()
  call slv%solve(1)
  print*, ""
  print '(3F6.3)', slv%state()

  call slv%dtor()
  !TODO: bcx, bcy, adv dtor

  contains
  subroutine init
    real, pointer :: tmp(:,:)
    tmp => slv%state() 
    print*, lbound(tmp)
    tmp = 0
    tmp(1:1,1:1) = 1

    tmp => slv%Cx() 
    tmp = .5

    tmp => slv%Cy() 
    tmp = .5
  end subroutine
end program
