!listing08
program test
  use solver_2D_m
  use mpdata_m
  use cyclic_m
  implicit none

  type(solver_2D_t) :: slv

  type(mpdata_t), target :: adv
  class(adv_t), pointer :: adv_p

  type(cyclic_t), target :: bcx, bcy
  class(bcd_t), pointer :: bcx_p, bcy_p

  integer :: nx = 3, ny = 2, nt = 4

  call adv%ctor
  call bcx%ctor
  call bcy%ctor
  adv_p => adv
  bcx_p => bcx
  bcy_p => bcy
  
  call slv%ctor(nx, ny, adv_p, bcx_p, bcy_p)

  call init
  call slv%solve(nt)
  call slv%dtor

  contains
  subroutine init
    real, dimension(:,:), pointer :: tmp
    tmp => slv%state() 
    tmp = 0
    tmp(0:2,0:3) = 1
  end subroutine
end program
!listing09
