module arr_m
  implicit none
  type :: arr_t
    real, dimension(:,:), pointer :: X
  end type
end module

module arakawa_c_m
  implicit none

  type :: half_t
  end type

  type(half_t) :: h

  interface operator (+)
    module procedure ph
  end interface

  interface operator (-)
    module procedure mh
  end interface

  ! TODO
  !interface operator (^)
  !  module procedure ext
  !end interface

  contains

  elemental function ph(i, h)
    integer, intent(in) :: i
    type(half_t), intent(in) :: h
    integer :: ph
    ph = i + 1
  end function

  elemental function mh(i, h)
    integer, intent(in) :: i 
    type(half_t), intent(in) :: h
    integer :: mh
    mh = i
  end function
end module

module adv_m
  use arr_m
  implicit none

  type, abstract :: adv_t
    integer :: n_steps, n_halos
    contains 
    procedure(op_2D_proto), deferred :: op_2D
  end type

  abstract interface
    subroutine op_2D_proto(this, psi, n, C, i, j, s)
      import :: adv_t
      import :: arr_t
      implicit none
      class(adv_t) :: this
      class(arr_t) :: psi(:), C(:)
      integer, intent(in) :: n, s
      integer, intent(in), dimension(:) :: i, j
    end subroutine
  end interface
end module

module bcd_m
  use arr_m
  implicit none

  type, abstract :: bcd_t
    contains
    procedure(fill_halos_proto), deferred :: fill_halos
  end type
 
  abstract interface 
    subroutine fill_halos_proto(this, psi)
      import :: arr_t
      import :: bcd_t
      implicit none
      class(bcd_t) :: this
      real, dimension(:,:), pointer :: psi
    end subroutine
  end interface
end module

module cyclic_m
  use bcd_m
  implicit none
  
  type, extends(bcd_t) :: cyclic_t
    contains
    procedure :: ctor => cyclic_ctor 
    procedure :: fill_halos => cyclic_fill_halos
  end type

  contains

  subroutine cyclic_ctor(this)
    implicit none
    class(cyclic_t) :: this
  end subroutine

  subroutine cyclic_fill_halos(this, psi)
    implicit none
    class(cyclic_t) :: this
    real, dimension(:,:), pointer :: psi
    !psi() = psi()
    !psi() = psi()
  end subroutine
end module

module mpdata_m
  use adv_m
  implicit none
  
  type, extends(adv_t) :: mpdata_t
    contains
    procedure :: ctor => mpdata_ctor
    procedure :: op_2D => mpdata_op_2D
  end type

  contains

  subroutine mpdata_ctor(this)
    implicit none
    class(mpdata_t) :: this
    this%n_steps = 2
    this%n_halos = 2
  end subroutine

  subroutine mpdata_op_2D(this, psi, n, C, i, j, s)
    use arakawa_c_m
    implicit none
    class(mpdata_t) :: this
    class(arr_t) :: psi(:), C(:)
    integer, intent(in) :: n, s
    integer, intent(in), dimension(:) :: i, j

    print*, "mpdata"
    psi(n+1)%X(i,j) = psi(n)%X(i,j) - ( &
      F(psi( n )%X( i,  j ), psi( n )%X( i+1,j ), C(0)%X( i + h, j )) - & 
      F(psi( n )%X( i-1,j ), psi( n )%X( i,  j ), C(0)%X( i - h, j ))   &   
    ) - ( &   
      F(psi( n )%X( i,j   ), psi( n )%X( i,j+1 ), C(1)%X( i,j + h )) - & 
      F(psi( n )%X( i,j-1 ), psi( n )%X( i,j   ), C(1)%X( i,j - h ))   &   
    )   

    contains 
    elemental function F(psi_l, psi_r, C)
      implicit none
      real :: F
      real, intent(in) :: psi_l, psi_r, C
    end function
  end subroutine
end module

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
    procedure :: ctor => solver_2D_ctor ! should be 'initial' but gfortran doesn't support it yet
    procedure :: solve => solver_2D_solve
    procedure :: state => solver_2D_state
    procedure :: dtor => solver_2D_dtor ! should be 'final' but gfortran doesn't support it yet
  end type 

  contains

  subroutine solver_2D_ctor(this, nx, ny, adv, bcx, bcy)
    use arakawa_c_m
    implicit none
    class(solver_2D_t) :: this
    integer, intent(in) :: nx, ny
    class(adv_t), pointer :: adv
    class(bcd_t), pointer :: bcx, bcy

    ! method body
    allocate(this%i(0:nx-1), this%j(0:ny-1))

    block
      integer :: i, hlo
      hlo = adv%n_halos
      
      allocate(this%psi(0:1))
      do i=0, 1
        allocate(this%psi(i)%X(- hlo - h: nx - 1 + hlo + h, - hlo : ny - 1 + hlo)) 
      end do

      allocate(this%C(0:1))
      do i=0, 1
        allocate(this%C(i)%X(- hlo : nx - 1 + hlo, - hlo - h : ny - 1 + hlo + h))   
      end do
    end block

    this%adv => adv
    this%bcx => bcx
    this%bcy => bcy
  end subroutine

  subroutine solver_2D_dtor(this)
    implicit none

    class(solver_2D_t) :: this

    integer :: i

    print*, 'destroying...'
    do i=0,1
      deallocate(this%psi(i)%X)   
    end do
    do i=0,1
      deallocate(this%C(i)%X)   
    end do
    deallocate(this%i, this%j, this%psi)
  end subroutine
  
  subroutine solver_2D_solve(this, nt)
    implicit none

    ! arguments
    class(solver_2D_t) :: this
    integer, intent(in) :: nt

    ! method body
    block
      integer :: t, s, n 
      n = 0
      s = 1
      do t = 1, nt
        !do s = 1, adv%n_steps
          print*, this%psi(0)%X(this%i,this%j)
          !this%bcx%fill_halos(this%psi(n))
          !this%bcy%fill_halos(this%psi(n))
          call this%adv%op_2D(this%psi, n, this%C, this%i, this%j, s)
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
