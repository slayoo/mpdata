module bcd_m
  use arr_m
  implicit none

  type, abstract :: bcd_t
    contains
    procedure(fill_halos_i), deferred :: fill_halos
  end type
 
  abstract interface 
    subroutine fill_halos_i(this, psi)
      import :: arr_t
      import :: bcd_t
      implicit none
      class(bcd_t) :: this
      real, dimension(:,:), pointer :: psi
    end subroutine
  end interface
end module
