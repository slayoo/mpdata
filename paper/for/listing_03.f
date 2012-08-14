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

    psi(n+1)%X(i,j) = psi(n)%X(i,j) - (    &
      F(                                   &
        psi(n)%X(i,  j), psi(n)%X(i+1,j),  &
        C(0)%X(i+h,j)                      &
      ) -                                  &
      F(                                   &
        psi(n)%X(i-1,j), psi(n)%X(i,  j),  &
        C(0)%X(i-h,j)                      &
      )                                    &
    ) - (                                  &
      F(                                   &
        psi(n)%X(i,j  ), psi(n)%X(i,j+1),  &
        C(1)%X(i,j+h)                      &
      ) -                                  &
      F(                                   &
        psi(n)%X(i,j-1), psi(n)%X(i,j  ),  &
        C(1)%X(i,j-h)                      &
      )                                    &
    )   

    contains 
    elemental function F(psi_l, psi_r, C)
      implicit none
      real :: F
      real, intent(in) :: psi_l, psi_r, C
    end function
  end subroutine
end module
