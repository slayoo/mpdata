module arr_m
  implicit none
  type :: arr_t
    real, dimension(:,:), pointer :: X
  end type
end module
