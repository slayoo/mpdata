program test_arrvec
  use arrvec_m
  class(arrvec_t), pointer :: psi
  integer, pointer, contiguous :: i(:), j(:)
  integer :: c, nx = 10

  allocate(i(nx), j(nx))
  i = (/ (c, c=0, nx-1) /)
  j = (/ (c, c=-2, nx-3) /)

  allocate(psi)
  call psi%ctor(2)
  call psi%init(0, i, j)
  call psi%init(1, i, j)

  print*, psi%at(0)%p%a
  print*, psi%at(0)%p%a(1,1)
  psi%at(0)%p%a(1,1) = 10
  print*, psi%at(0)%p%a(1,1)
  print*, psi%at(-2)%p%a(1,1)

  call psi%dtor
  deallocate(psi)
  deallocate(i, j)
end

