!listing00
! code licensed under the terms of GNU GPL v3
! copyright holder: University of Warsaw
!listing01
module real_m
  implicit none
  integer, parameter :: real_t = kind(0.d0) 
end module
!listing02
module arrvec_m
  use real_m
  implicit none

  type :: arr_t
    real(real_t), pointer, contiguous :: a(:,:)
  end type

  type :: arrptr_t
    class(arr_t), pointer :: p
  end type

  type :: arrvec_t
    class(arrptr_t), pointer :: at(:)
    logical, pointer, contiguous :: inited(:)
    contains
    procedure :: ctor => arrvec_ctor
    procedure :: init => arrvec_init
    procedure :: dtor => arrvec_dtor
  end type

  contains

  subroutine arrvec_ctor(O, n)
    class(arrvec_t) :: O
    integer, intent(in) :: n

    allocate(O%at(-n:n-1))
    allocate(O%inited(0:n-1))
    O%inited = .false.
  end subroutine

  subroutine arrvec_init(O, n,    &
    i_min, i_max, j_min, j_max       &
  )
    class(arrvec_t) :: O
    integer, intent(in) :: n,        &
      i_min, i_max, j_min, j_max

    allocate(O%at(n)%p)
    allocate(O%at(n)%p%a(         &
      i_min : i_max, j_min : j_max   &
    ))
    O%inited(n) = .true.
    O%at(n - size(O%inited))%p &
      => O%at(n)%p
  end subroutine

  subroutine arrvec_dtor(O)
    class(arrvec_t) :: O
    integer :: i

    do i = 0, size(O%inited) - 1
      if (O%inited(i)) then
        deallocate(O%at(i)%p%a)
        deallocate(O%at(i)%p)
      end if
    end do
    deallocate(O%at)
  end subroutine
end module
!listing03
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

  contains

  elemental function ph(i, h) result (return)
    integer, intent(in) :: i
    type(half_t), intent(in) :: h
    integer :: return
    return = i + 1
  end function

  elemental function mh(i, h) result (return)
    integer, intent(in) :: i 
    type(half_t), intent(in) :: h
    integer :: return
    return = i
  end function
end module
!listing04
module adv_m
  use arrvec_m
  implicit none

  type, abstract :: adv_t
    integer :: n_steps, n_halos
    contains 
    procedure(op_2D_i), deferred :: op_2D
    procedure(adv_init_i), deferred :: init
  end type

  abstract interface
    subroutine op_2D_i(O, psi, n, C, step)
      import :: adv_t
      import :: arrvec_t
      class(adv_t) :: O
      class(arrvec_t), pointer :: psi, C
      integer, intent(in) :: n, step
    end subroutine

    subroutine adv_init_i(O, i, j)
      import :: adv_t
      class(adv_t) :: O
      integer, intent(in), contiguous, pointer :: i(:), j(:) 
    end subroutine
  end interface
end module
!listing05
module bcd_m
  use arrvec_m
  implicit none

  type, abstract :: bcd_t
    contains
    procedure(bcd_fill_halos_i), deferred :: fill_halos
    procedure(bcd_init_i), deferred :: init
  end type
 
  abstract interface 
    subroutine bcd_fill_halos_i(O, psi)
      import :: bcd_t, real_t
      class(bcd_t) :: O
      real(real_t), pointer, contiguous :: psi(:,:) 
    end subroutine
    subroutine bcd_init_i(O, d, i, j, hlo)
      import :: bcd_t, real_t
      class(bcd_t) :: O
      integer :: d, hlo
      integer, pointer, contiguous :: i(:), j(:)
    end subroutine
  end interface
end module
!listing06
module solver_2D_m
  use arrvec_m
  use adv_m
  use bcd_m
  implicit none

  type :: solver_2D_t
    class(arrvec_t), pointer :: psi, C
    class(adv_t), pointer :: adv
    class(bcd_t), pointer :: bcx, bcy
    integer, pointer, contiguous :: i(:), j(:) 
    integer :: n
    contains
    procedure :: ctor  => solver_2D_ctor 
    procedure :: solve => solver_2D_solve
    procedure :: state => solver_2D_state
    procedure :: Cx => solver_2D_Cx
    procedure :: Cy => solver_2D_Cy
    procedure :: dtor  => solver_2D_dtor 
  end type 

  contains

  subroutine solver_2D_ctor(O, nx,ny,adv,bcx,bcy)
    use arakawa_c_m
    class(solver_2D_t) :: O
    integer, intent(in) :: nx, ny
    class(adv_t), intent(in), pointer :: adv
    class(bcd_t), intent(in), pointer :: bcx, bcy

    allocate(O%i(0:nx-1), O%j(0:ny-1))
    block
      integer :: c
      O%i = (/ (c, c=0, nx-1) /)
      O%j = (/ (c, c=0, ny-1) /)
    end block

    call bcx%init(0, O%i, O%j, adv%n_halos)
    call bcy%init(1, O%j, O%i, adv%n_halos)
    call adv%init(O%i, O%j)

    block
      integer :: i, hlo
      hlo = adv%n_halos

      allocate(O%psi)
      call O%psi%ctor(2)
      
      do i=0, 1
        call O%psi%init(i, &
          O%i(0) -hlo, O%i(size(O%i)-1) + hlo, &
          O%j(0) -hlo, O%j(size(O%j)-1) + hlo  &
        ) 
      end do

      allocate(O%C)
      call O%C%ctor(2)
      call O%C%init(0,                                     &
        O%i(0) -hlo - h, O%i(size(O%i)-1) + hlo + h, &
        O%j(0) -hlo    , O%j(size(O%j)-1) + hlo      &   
      )
      call O%C%init(1,                                     &
        O%i(0) -hlo    , O%i(size(O%i)-1) + hlo,     &   
        O%j(0) -hlo - h, O%j(size(O%j)-1) + hlo + h  &
      )
    end block

    O%n = 0

    O%adv => adv
    O%bcx => bcx
    O%bcy => bcy
  end subroutine

  subroutine solver_2D_dtor(O)
    class(solver_2D_t) :: O
    call O%psi%dtor()
    deallocate(O%i, O%j, O%psi)
  end subroutine
  
  subroutine solver_2D_solve(O, nt)
    class(solver_2D_t) :: O
    integer, intent(in) :: nt

    block
      integer :: t, s 
      do t = 0, nt-1
        do s = 0, O%adv%n_steps-1
          call O%bcx%fill_halos(O%psi%at(O%n)%p%a)
          call O%bcy%fill_halos(O%psi%at(O%n)%p%a)
          call O%adv%op_2D(O%psi, O%n, O%C, s)
          O%n = mod(O%n + 1 + 2, 2) - 2
        end do
      end do
    end block
  end subroutine

  function solver_2D_state(O) result (return)
    class(solver_2D_t) :: O
    real(real_t), pointer :: return(:,:)
    return => O%psi%at(O%n)%p%a( &
      O%i(0) : O%i(size(O%i)-1), &
      O%j(0) : O%j(size(O%j)-1)  &
    )
  end function

  function solver_2D_Cx(O) result (return)
    class(solver_2D_t) :: O
    real(real_t), pointer :: return(:,:)
    return => O%C%at(0)%p%a
  end function

  function solver_2D_Cy(O) result (return)
    class(solver_2D_t) :: O
    real(real_t), pointer :: return(:,:)
    return => O%C%at(1)%p%a
  end function
end module
!listing
module pi_m
  use real_m
  contains
  function pi(d, arr, i, j) result(return)
    integer :: d
    real(real_t), pointer :: arr(:,:), return(:,:)
    integer :: i(:), j(:)
    select case (d) 
      case (0) 
        return => arr(        &   
          i(1) : i(size(i)),  &
          j(1) : j(size(j))   &   
        )   
      case (1) 
        return => arr(        &   
          j(1) : j(size(j)),  &
          i(1) : i(size(i))   &   
        )   
    end select
  end function
end module
!listing07
module cyclic_m
  use bcd_m
  use adv_m
  use pi_m
  implicit none
  
  type, extends(bcd_t) :: cyclic_t
    integer :: d
    integer, pointer :: left_halo(:), rght_halo(:)
    integer, pointer :: left_edge(:), rght_edge(:)
    integer, pointer, contiguous :: j(:)
    contains
    procedure :: init => cyclic_init
    procedure :: dtor => cyclic_dtor 
    procedure :: fill_halos => cyclic_fill_halos
  end type

  contains

  subroutine cyclic_init(O, d, i, j, hlo)
    class(cyclic_t) :: O
    integer :: d, hlo
    integer, pointer, contiguous :: i(:), j(:) 

    O%d = d

    allocate(O%left_halo(hlo))
    allocate(O%rght_halo(hlo))
    allocate(O%left_edge(hlo))
    allocate(O%rght_edge(hlo))
    allocate(O%j(size(j) + 2*hlo))

    block
      integer :: c
      O%left_halo = (/ (c, c = i(0) - hlo, i(0) - 1) /) 
      O%rght_halo = (/ (c, c = i(size(i)-1) + 1, i(size(i)-1) + hlo) /) 
      O%left_edge = (/ (c, c = i(0),  i(0) + hlo) /)
      O%rght_edge = (/ (c, c = i(size(i) - hlo), i(size(i)-1)) /)
      O%j         = (/ (c, c = j(0) - hlo, j(size(j)-1) + hlo) /)
    end block
  end subroutine

  subroutine cyclic_dtor(O)
    class(cyclic_t) :: O
    deallocate(O%left_halo)
    deallocate(O%rght_halo)
    deallocate(O%left_edge)
    deallocate(O%rght_edge)
  end subroutine

  subroutine cyclic_fill_halos(O, psi)
    class(cyclic_t) :: O
    real(real_t), pointer :: psi(:,:), tmp(:,:)
    tmp => pi(O%d, psi, O%left_halo, O%j) 
    tmp = pi(O%d, psi, O%rght_edge, O%j)
    tmp => pi(O%d, psi, O%rght_halo, O%j) 
    tmp = pi(O%d, psi, O%left_edge, O%j)
  end subroutine
end module
!listing08
module donorcell_1D_m
  use real_m
  use arakawa_c_m
  use pi_m
  implicit none
  contains 

  elemental function F(psi_l, psi_r, C) result (return)
    real(real_t) :: return
    real(real_t), intent(in) :: psi_l, psi_r, C
    return = (                            &
      (C + abs(C)) * psi_l +              &
      (C - abs(C)) * psi_r                &
    ) / 2
  end function
 
  function donorcell(d, psi, C, i, j) result (return)
    integer :: d
    integer, pointer, intent(in), contiguous :: i(:), j(:) 
    real(real_t) :: return(size(i), size(j))
    real(real_t), pointer, intent(in), contiguous :: psi(:,:), C(:,:)
    return = (                                 &
      F(pi(d, psi, i,   j), pi(d, psi, i+1, j), pi(d, C, i+h, j)) - &
      F(pi(d, psi, i-1, j), pi(d, psi, i,   j), pi(d, C, i-h, j))   &
    )
  end function
end module
!listing09
module donorcell_2D_m
  use donorcell_1D_m
  use arrvec_m
  implicit none
  contains 
  subroutine donorcell_2D(psi, n, C, i, j)  
    class(arrvec_t), pointer :: psi, C
    integer, intent(in) :: n
    integer, pointer, intent(in), contiguous :: i(:), j(:) 
    
    psi%at( n+1 )%p%a( i, j ) = psi%at( n )%p%a( i, j ) &
      - donorcell(0, psi%at(n)%p%a, C%at(0)%p%a, i, j)  &
      - donorcell(1, psi%at(n)%p%a, C%at(1)%p%a, j, i)      
  end subroutine
end module
!listing10
module antidiff_2D_m
  use arrvec_m
  use arakawa_c_m
  use pi_m
  implicit none

  contains 

  function frac(nom, den) result (return)
    real(real_t), intent(in) :: nom(:,:), den(:,:)
    real(real_t) :: return(size(nom,1), size(nom,2))
    where (den > 0)
      return = nom / den
    elsewhere
      return = 0
    end where
  end function

  function A(d, psi, i, j) result (return)
    integer :: d
    real(real_t), pointer, intent(in), contiguous :: psi(:,:)
    integer, pointer, contiguous :: i(:), j(:) 
    real(real_t) :: return(size(i), size(j))
    return = frac(                           &
      pi(d, psi, i+1, j) - pi(d, psi, i, j), &
      pi(d, psi, i+1, j) + pi(d, psi, i, j)  &
    )  
  end function

  function B(d, psi, i, j) result (return)
    integer :: d
    real(real_t), pointer, intent(in), contiguous :: psi(:,:) 
    integer, pointer, contiguous :: i(:), j(:) 
    real(real_t) :: return(size(i), size(j))
    return = frac(   &
      pi(d, psi, i+1, j+1) + pi(d, psi, i,   j+1)  &
    - pi(d, psi, i+1, j-1) - pi(d, psi, i,   j-1), &
      pi(d, psi, i+1, j+1) + pi(d, psi, i,   j+1)  &
    + pi(d, psi, i+1, j-1) + pi(d, psi, i,   j-1)  &
    ) / 2
  end function

  function antidiff_2D(d, psi, i, j, C) result (return)
    integer :: d
    integer, pointer, contiguous :: i(:), j(:) 
    real(real_t) :: return(size(i), size(j))
    real(real_t), pointer, intent(in), contiguous :: psi(:,:) 
    class(arrvec_t), pointer :: C
 
    return =                                  &
      abs(pi(d, C%at(d)%p%a, i+h, j))         &
      * (1 - abs(pi(d, C%at(d)%p%a, i+h, j))) &
      * A(d, psi, i, j)                       &
      - pi(d, C%at(d)%p%a, i+h, j)            &
      * (                                     &
        pi(d, C%at(d-1)%p%a, i+1, j+h) +      &
        pi(d, C%at(d-1)%p%a, i,   j+h) +      &
        pi(d, C%at(d-1)%p%a, i+1, j-h) +      &
        pi(d, C%at(d-1)%p%a, i,   j-h)        &
      ) / 4                                   &
      * B(d, psi, i, j)
  end function
end module
!listing08
module mpdata_m
  use adv_m
  use donorcell_2D_m
  use antidiff_2D_m
  implicit none
  
  type, extends(adv_t) :: mpdata_t
    integer :: n_iters, n_tmp
    class(arrvec_t), pointer :: tmp(:) 
    integer, pointer, contiguous :: i(:), j(:)
    integer, pointer, contiguous :: im(:), jm(:) 
    contains
    procedure :: ctor => mpdata_ctor
    procedure :: init => mpdata_init
    procedure :: op_2D => mpdata_op_2D
  end type

  contains

  subroutine mpdata_ctor(O, n_iters)
    class(mpdata_t) :: O
    integer, intent(in) :: n_iters
    O%n_iters = n_iters
    O%n_steps = n_iters
    O%n_halos = n_iters
  end subroutine 

  subroutine mpdata_init(O, i, j)
    class(mpdata_t) :: O
    integer, intent(in), pointer, contiguous :: i(:), j(:)

    O%i => i
    O%j => j

    allocate(O%im(0:size(i)), O%jm(0:size(j)))  
    block 
      integer :: c
      O%im = (/ (c, c=i(0)-1, j(size(i)-1)) /)
      O%jm = (/ (c, c=j(0)-1, j(size(j)-1)) /)
    end block
    
    block
      integer :: c, hlo, nc
      O%n_tmp = 1
      if (O%n_iters > 2) O%n_tmp = 2
      allocate(O%tmp(0:O%n_tmp)) 
      hlo = O%n_halos
      do c=0, O%n_tmp - 1
        call O%tmp(c)%ctor(2)
        call O%tmp(c)%init(0,                    &   
          i(0) -hlo - h, i(size(i)-1) + hlo + h, &
          j(0) -hlo    , j(size(j)-1) + hlo      &   
        )   
        call O%tmp(c)%init(1,                    &   
          i(0) -hlo    , i(size(i)-1) + hlo,     &   
          j(0) -hlo - h, j(size(j)-1) + hlo + h  &
        )  
      end do
    end block
  end subroutine

  subroutine mpdata_dtor(O)
    class(mpdata_t) :: O
    integer :: c
    do c=0, O%n_tmp-1 
      call O%tmp(c)%dtor()
    end do
    deallocate(O%im, O%jm, O%tmp)
  end subroutine

  subroutine mpdata_op_2D(O, psi, n, C, step)
    class(mpdata_t) :: O
    class(arrvec_t), pointer :: psi, C
    integer, intent(in) :: n, step
    if (step == 0) then
      call donorcell_2D(psi, n, C, O%i, O%j)
    else
      block
        class(arrvec_t), pointer :: C_corr, C_unco
        if (step == 1) then
          C_unco => C
          C_corr => O%tmp(0)
        else if (mod(step, 2) == 1) then
          C_unco => O%tmp(1)
          C_corr => O%tmp(0)
        else
          C_unco => O%tmp(0)
          C_corr => O%tmp(1)
        end if

        !TODO: hlo!
        C_corr%at(0)%p%a(O%im+h, O%j) = &
          antidiff_2D(0, psi%at(n)%p%a, O%im, O%j, C_unco)
        C_corr%at(1)%p%a(O%i, O%jm+h) = &
          antidiff_2D(1, psi%at(n)%p%a, O%jm, O%i, C_unco)

        call donorcell_2D(psi, n, C_corr, O%i, O%j)
      end block
    end if
  end subroutine
end module
!listing09
