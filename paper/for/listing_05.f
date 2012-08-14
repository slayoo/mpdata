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
