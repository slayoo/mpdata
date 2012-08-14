module adv_m
  use arr_m
  implicit none

  type, abstract :: adv_t
    integer :: n_steps, n_halos
    contains 
    procedure(op_2D_i), deferred :: op_2D
  end type

  abstract interface
    subroutine op_2D_i(this, psi, n, C, i, j, s)
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
