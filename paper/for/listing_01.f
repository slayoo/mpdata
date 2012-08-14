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
