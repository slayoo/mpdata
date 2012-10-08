!listing00
module arrvec_m
  implicit none

  type :: arr_t
    real, pointer :: a(:,:)
  end type

  type :: arrptr_t
    class(arr_t), pointer :: p
  end type

  type :: arrvec_t
    class(arrptr_t), pointer :: at(:)
    logical, pointer :: inited(:)
    contains
    procedure :: ctor => arrvec_ctor
    procedure :: init => arrvec_init
    procedure :: dtor => arrvec_dtor
  end type

  contains

  subroutine arrvec_ctor(this, n)
    class(arrvec_t) :: this
    integer, intent(in) :: n

    allocate(this%at(-n:n-1))
    allocate(this%inited(0:n-1))
    this%inited = .false.
  end subroutine

  subroutine arrvec_init(this, n, i_min, i_max, j_min, j_max)
    class(arrvec_t) :: this
    integer, intent(in) :: n, i_min, i_max, j_min, j_max

    allocate(this%at(n)%p)
    allocate(this%at(n)%p%a(i_min : i_max, j_min : j_max))
    this%inited(n) = .true.
    this%at(n - size(this%inited))%p => this%at(n)%p
  end subroutine

  subroutine arrvec_dtor(this)
    class(arrvec_t) :: this
    integer :: i

    do i = 0, size(this%inited) - 1
      if (this%inited(i)) then
        deallocate(this%at(i)%p%a)
        deallocate(this%at(i)%p)
      end if
    end do
    deallocate(this%at)
  end subroutine
end module
!listing01
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
!listing02
module adv_m
  use arrvec_m
  implicit none

  type, abstract :: adv_t
    integer :: n_steps, n_halos
    contains 
    procedure(op_2D_i), deferred :: op_2D
  end type

  abstract interface
    subroutine op_2D_i(this, psi, n, C, i, j, s)
      import :: adv_t
      import :: arrvec_t
      class(adv_t) :: this
      class(arrvec_t), pointer :: psi, C
      integer, intent(in) :: n, s
      integer, intent(in) :: i(:), j(:)
    end subroutine
  end interface
end module
!listing03
module bcd_m
  use arrvec_m
  implicit none

  type, abstract :: bcd_t
    contains
    procedure(fill_halos_i), deferred :: fill_halos
  end type
 
  abstract interface 
    subroutine fill_halos_i(this, psi)
      import :: arrvec_t
      import :: bcd_t
      class(bcd_t) :: this
      real, pointer :: psi(:,:)
    end subroutine
  end interface
end module
!listing04
module solver_2D_m
  use arrvec_m
  use adv_m
  use bcd_m
  implicit none

  type :: solver_2D_t
    class(arrvec_t), pointer :: psi, C
    class(adv_t), pointer :: adv
    class(bcd_t), pointer :: bcx, bcy
    integer, pointer :: i(:), j(:)
    integer :: n
    contains
    procedure :: ctor  => solver_2D_ctor 
    procedure :: solve => solver_2D_solve
    procedure :: state => solver_2D_state
    procedure :: dtor  => solver_2D_dtor 
  end type 

  contains

  subroutine solver_2D_ctor(this, nx,ny,adv,bcx,bcy)
    use arakawa_c_m
    class(solver_2D_t) :: this
    integer, intent(in) :: nx, ny
    class(adv_t), pointer :: adv
    class(bcd_t), pointer :: bcx, bcy

    allocate(this%i(0:nx-1), this%j(0:ny-1))

    block
      integer :: i, hlo
      hlo = adv%n_halos

      allocate(this%psi)
      call this%psi%ctor(2)
      
      do i=0, 1
        call this%psi%init(i,         &
          -hlo - h, nx - 1 + hlo + h, &
          -hlo    , ny - 1 + hlo      &
        ) 
      end do

      allocate(this%C)
      call this%C%ctor(2)
      do i=0, 1 !TODO: different dims!!!
        call this%C%init(i,           &
          -hlo    , nx - 1 + hlo,     &
          -hlo - h, ny - 1 + hlo + h  &
        )
      end do
    end block

    this%n = 0

    this%adv => adv
    this%bcx => bcx
    this%bcy => bcy
  end subroutine

  subroutine solver_2D_dtor(this)
    class(solver_2D_t) :: this
    call this%psi%dtor
    deallocate(this%i, this%j, this%psi)
  end subroutine
  
  subroutine solver_2D_solve(this, nt)
    class(solver_2D_t) :: this
    integer, intent(in) :: nt

    block
      integer :: t, s 
      s = 1
      do t = 1, nt
        !do s = 1, adv%n_steps
          !this%bcx%fill_halos(this%psi(n))
          !this%bcy%fill_halos(this%psi(n))
          call this%adv%op_2D( &
            this%psi, this%n, this%C, this%i, this%j, s &
          )
          this%n = mod(this%n + 1 + 2, 2) - 2
        !end do
      end do
    end block
  end subroutine

  function solver_2D_state(this) result (state)
    class(solver_2D_t) :: this
    real, pointer :: state(:,:)
    state => this%psi%at(this%n)%p%a( &
      this%i(0) : this%i(size(this%i)-1), &
      this%j(0) : this%j(size(this%j)-1)  &
    )
  end function
end module
!listing05
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
    class(cyclic_t) :: this
  end subroutine

  subroutine cyclic_fill_halos(this, psi)
    class(cyclic_t) :: this
    real, pointer :: psi(:,:)
    !psi() = psi()
    !psi() = psi()
  end subroutine
end module
!listing06
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
    class(mpdata_t) :: this
    this%n_steps = 2
    this%n_halos = 2
  end subroutine

  subroutine mpdata_op_2D(this, psi, n, C, i, j, s)
    use arakawa_c_m
    class(mpdata_t) :: this
    class(arrvec_t), pointer :: psi, C
    integer, intent(in) :: n, s
    integer, intent(in) :: i(:), j(:)

    psi%at(n+1)%p%a(i,j) = psi%at(n)%p%a(i,j) - (    &
      F(                                             &
        psi%at( n )%p%a( i,   j ),  &
        psi%at( n )%p%a( i+1, j ),  &
        C%at( 0 )%p%a( i + h, j )   &
      ) -                                            &
      F(                                             &
        psi%at( n )%p%a( i-1, j ),  &
        psi%at( n )%p%a( i,   j ),  &
        C%at( 0 )%p%a( i-h, j )     &
      )                                              &
    ) - (                                            &
      F(                                             &
        psi%at(n)%p%a(i,j  ), psi%at(n)%p%a(i,j+1),  &
        C%at(1)%p%a(i,j+h)                           &
      ) -                                            &
      F(                                             &
        psi%at(n)%p%a(i,j-1), psi%at(n)%p%a(i,j  ),  &
        C%at(1)%p%a(i,j-h)                           &
      )                                              &
    )   

    contains 
    elemental function F(psi_l, psi_r, C)
      real :: F
      real, intent(in) :: psi_l, psi_r, C
    end function
  end subroutine
end module
!listing07
