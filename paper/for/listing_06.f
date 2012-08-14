module solver_2D_m
  use arr_m
  use adv_m
  use bcd_m
  implicit none

  type :: solver_2D_t
    class(arr_t), allocatable :: psi(:), C(:)
    class(adv_t), pointer :: adv
    class(bcd_t), pointer :: bcx, bcy
    integer,  dimension(:), allocatable :: i, j
    contains
    procedure :: ctor  => solver_2D_ctor 
    procedure :: solve => solver_2D_solve
    procedure :: state => solver_2D_state
    procedure :: dtor  => solver_2D_dtor 
  end type 

  contains

  subroutine solver_2D_ctor(this, nx,ny,adv,bcx,bcy)
    use arakawa_c_m
    implicit none
    class(solver_2D_t) :: this
    integer, intent(in) :: nx, ny
    class(adv_t), pointer :: adv
    class(bcd_t), pointer :: bcx, bcy

    allocate(this%i(0:nx-1), this%j(0:ny-1))

    block
      integer :: i, hlo
      hlo = adv%n_halos
      
      allocate(this%psi(0:1))
      do i=0, 1
        allocate(this%psi(i)%X(        &
          -hlo - h : nx - 1 + hlo + h, &
          -hlo     : ny - 1 + hlo      &
        )) 
      end do

      allocate(this%C(0:1))
      do i=0, 1
        allocate(this%C(i)%X(          &
          -hlo     : nx - 1 + hlo,     &
          -hlo - h : ny - 1 + hlo + h  &
        ))   
      end do
    end block

    this%adv => adv
    this%bcx => bcx
    this%bcy => bcy
  end subroutine

  subroutine solver_2D_dtor(this)
    implicit none
    class(solver_2D_t) :: this

    block
      integer :: i
      do i=0,1
        deallocate(this%psi(i)%X)   
      end do
      do i=0,1
        deallocate(this%C(i)%X)   
      end do
    end block
    deallocate(this%i, this%j, this%psi)
  end subroutine
  
  subroutine solver_2D_solve(this, nt)
    implicit none
    class(solver_2D_t) :: this
    integer, intent(in) :: nt

    block
      integer :: t, s, n 
      n = 0
      s = 1
      do t = 1, nt
        !do s = 1, adv%n_steps
          !this%bcx%fill_halos(this%psi(n))
          !this%bcy%fill_halos(this%psi(n))
          call this%adv%op_2D( &
            this%psi, n, this%C, this%i, this%j, s &
          )
        !end do
      end do
    end block
  end subroutine

  function solver_2D_state(this) result (state)
    class(solver_2D_t) :: this
    real, dimension(:,:), pointer :: state
    state => this%psi(0)%X( &
      this%i(0) : this%i(size(this%i)-1), &
      this%j(0) : this%j(size(this%j)-1)  &
    )
  end function
end module
