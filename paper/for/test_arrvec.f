program test_arrvec
  use arrvec_m
  class(arrvec_t), pointer :: psi

  allocate(psi)
  call psi%ctor(2)
  call psi%init(0, 0, 3, 0, 4)

  print*, psi%at(0)%p%a
  print*, psi%at(0)%p%a(1,1)
  psi%at(0)%p%a(1,1) = 10
  print*, psi%at(0)%p%a(1,1)

  call psi%dtor
  deallocate(psi)
end

