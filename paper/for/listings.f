!listing00
! code licensed under the terms of GNU GPL v3
! copyright holder: University of Warsaw
!listing01
module real_m
  implicit none
  integer,parameter:: real_t = kind(0.d0) 
end module
!listing02
module arrvec_m
  use real_m
  implicit none

  type:: arr_t
    real(real_t), pointer, contiguous :: a(:,:)
  end type

  type:: arrptr_t
    class(arr_t), pointer :: p
  end type

  type:: arrvec_t
    class(arrptr_t), pointer :: at(:)
    logical, pointer, contiguous :: inited(:)
    contains
    procedure :: ctor => arrvec_ctor
    procedure :: init => arrvec_init
    procedure :: dtor => arrvec_dtor
  end type

  contains

  subroutine arrvec_ctor(this, n)
    class(arrvec_t):: this
    integer, intent(in) :: n
    allocate(this%at( -n : n-1 ))
    allocate(this%inited( 0 : n-1 ))
    this%inited = .false.
  end subroutine

  subroutine arrvec_init(this, n, i, j)
    class(arrvec_t) :: this
    integer, intent(in) :: n
    integer, intent(in), dimension(:) :: i, j
    allocate(this%at(n)%p)
    allocate(this%at(n)%p%a(                  &
      i(1) : i(size(i)),            &
      j(1) : j(size(j))             &
    ))
    this%inited(n) = .true.
    this%at(n - size(this%inited))%p          &
      => this%at(n)%p
  end subroutine

  subroutine arrvec_dtor(this)
    class(arrvec_t):: this
    integer:: i
    do i = 0, size(this%inited) - 1
      if (this%inited(i)) then
        deallocate(this%at(i)%p%a)
        deallocate(this%at(i)%p)
      end if
    end do
    deallocate(this%at)
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
    return = i 
  end function

  elemental function mh(i, h) result (return)
    integer, intent(in) :: i 
    type(half_t), intent(in) :: h
    integer :: return
    return = i - 1
  end function
end module
!
module halo_m
  use arakawa_c_m
  implicit none

  interface operator (//)
    module procedure pmn
    module procedure pmh
  end interface 

  contains

  function pmn(r, n) result (return)
    integer, intent(in), dimension(:) :: r
    integer, intent(in) :: n
    integer, dimension(size(r)+2) :: return
    
    integer :: c
    return = (/ (c, c=r(1) - n, r(size(r)) + n) /)
  end function

  function pmh(r, h) result (return)
    integer, intent(in), dimension(:) :: r
    type(half_t), intent(in) :: h
    integer, dimension(size(r)+1) :: return
    
    integer :: c
    return = (/ (c, c=r(1) - h, r(size(r)) + h) /)
  end function
end module
!listing04
module pi_m
  use real_m
  contains
  function pi(d, arr, i, j) result(return)
    integer:: d
    real(real_t),pointer:: arr(:,:), return(:,:)
    integer:: i(:), j(:)
    select case (d) 
      case (0) 
        return => arr(                        &   
          i(1) : i(size(i)),                  &
          j(1) : j(size(j))                   &   
        )   
      case (1) 
        return => arr(                        &   
          j(1) : j(size(j)),                  &
          i(1) : i(size(i))                   &   
        )   
    end select
  end function
end module
!listing05
module bcd_m
  use arrvec_m
  implicit none

  type,abstract:: bcd_t
    contains
    procedure(bcd_fill_halos), deferred :: fill_halos
    procedure(bcd_init), deferred :: init
  end type
 
  abstract interface 
    subroutine bcd_fill_halos(this, a, j)
      import :: bcd_t, real_t
      class(bcd_t ):: this
      real(real_t), pointer, contiguous :: a(:,:) 
      integer, dimension(:) :: j
    end subroutine
    subroutine bcd_init(this, d, n, hlo)
      import :: bcd_t
      class(bcd_t) :: this
      integer :: d, n, hlo
    end subroutine
  end interface
end module
!listing06
module solver_2D_m
  use arrvec_m
  use bcd_m
  use arakawa_c_m
  use halo_m
  implicit none

  type,abstract:: solver_2D_t
    class(arrvec_t),pointer:: psi, C
    integer:: n, hlo
    integer,pointer,contiguous:: i(:), j(:) 
    class(bcd_t),pointer:: bcx, bcy
    contains
    procedure:: solve   => solver_2D_solve
    procedure:: state   => solver_2D_state
    procedure:: courant => solver_2D_courant
    procedure:: cycle   => solver_2D_cycle
    procedure(solver_2D_advop),deferred:: advop
  end type 

  abstract interface
    subroutine solver_2D_advop(this)
      import solver_2D_t
      class(solver_2D_t):: this
    end subroutine
  end interface

  contains

  subroutine solver_2D_ctor(this, bcx, bcy, nx, ny, hlo)
    use arakawa_c_m
    use halo_m
    class(solver_2D_t):: this
    class(bcd_t),intent(in),pointer:: bcx,bcy
    integer,intent(in):: nx, ny, hlo

    this%n = 0
    this%hlo = hlo
    this%bcx => bcx
    this%bcy => bcy

    allocate(this%i(0:nx-1), this%j(0:ny-1))
    block
      integer:: c
      this%i = (/ (c, c=0, nx-1) /)
      this%j = (/ (c, c=0, ny-1) /)
    end block

    call bcx%init(0, nx, hlo)
    call bcy%init(1, ny, hlo)

    allocate(this%psi)
    call this%psi%ctor(2)
    block
      integer:: n
      do n=0, 1
        call this%psi%init(n, this%i // hlo, this%j // hlo)
      end do
    end block

    allocate(this%C)
    call this%C%ctor(2)
    call this%C%init(0, this%i // h, this%j // hlo)
    call this%C%init(1, this%i // hlo, this%j // h)
  end subroutine

  subroutine solver_2D_dtor(this)
    class(solver_2D_t):: this
    call this%psi%dtor()
    call this%C%dtor()
    deallocate( &
      this%i,   &
      this%j,   &
      this%psi, &
      this%C    &
    )
  end subroutine
  
  function solver_2D_state(this) result (return)
    class(solver_2D_t):: this
    real(real_t),pointer:: return(:,:)
    return => this%psi%at(this%n)%p%a( &
      this%i(0) : this%i(size(this%i)-1), &
      this%j(0) : this%j(size(this%j)-1)  &
    )
  end function

  function solver_2D_courant(this, d) result (return)
    class(solver_2D_t):: this
    integer :: d
    real(real_t),pointer :: return(:,:)
    return => this%C%at(d)%p%a(               & 
      this%i(0)-h : this%i(size(this%i)-1)+h, &
      this%j(0)-h : this%j(size(this%j)-1)+h  &
    )
  end function

  subroutine solver_2D_cycle(this)
    class(solver_2D_t):: this
    this%n = mod(this%n + 1 + 2, 2) - 2
  end subroutine

  subroutine solver_2D_solve(this, nt)
    class(solver_2D_t) :: this
    integer, intent(in) :: nt
    integer :: t

    do t = 0, nt-1 
      call this%bcx%fill_halos(this%psi%at(this%n)%p%a, this%j // this%hlo)
      call this%bcy%fill_halos(this%psi%at(this%n)%p%a, this%i // this%hlo)
      call this%advop()
      call this%cycle()
    end do
  end subroutine
end module
!listing07
module cyclic_m
  use bcd_m
  use pi_m
  implicit none
  
  type, extends(bcd_t) :: cyclic_t
    integer :: d
    integer, pointer, contiguous :: left_halo(:), rght_halo(:), left_edge(:), rght_edge(:) 
    contains
    procedure:: init => cyclic_init
    procedure:: dtor => cyclic_dtor 
    procedure:: fill_halos => cyclic_fill_halos
  end type

  contains

  subroutine cyclic_init(this, d, n, hlo)
    class(cyclic_t):: this
    integer:: d, n, hlo

    this%d = d
    allocate(this%left_halo(hlo))
    allocate(this%rght_halo(hlo))
    allocate(this%left_edge(hlo))
    allocate(this%rght_edge(hlo))

    block
      integer :: c
      this%left_halo = (/(c, c=-hlo, -1)/) 
      this%rght_halo = (/(c, c=n, n-1+hlo)/) 
      this%left_edge = (/(c, c=0, hlo-1)/)
      this%rght_edge = (/(c, c=n-hlo, n-1)/)
    end block
  end subroutine

  subroutine cyclic_dtor(this)
    class(cyclic_t):: this
    deallocate(this%left_halo)
    deallocate(this%rght_halo)
    deallocate(this%left_edge)
    deallocate(this%rght_edge)
  end subroutine

  subroutine cyclic_fill_halos(this, a, j)
    class(cyclic_t):: this
    real(real_t), pointer :: a(:,:), tmp(:,:)
    integer, dimension(:) :: j
    tmp => pi(this%d, a, this%left_halo, j) 
    tmp =  pi(this%d, a, this%rght_edge, j)
    tmp => pi(this%d, a, this%rght_halo, j) 
    tmp =  pi(this%d, a, this%left_edge, j)
  end subroutine
end module
!listing08
module donorcell_m
  use real_m
  use arakawa_c_m
  use pi_m
  use arrvec_m
  implicit none
  contains 
!listing09
  elemental function F(psi_l, psi_r, C) result (return)
    real(real_t):: return
    real(real_t),intent(in):: psi_l, psi_r, C
    return = (                                &
      (C + abs(C)) * psi_l +                  &
      (C - abs(C)) * psi_r                    &
    ) / 2
  end function
!listing10
  function donorcell(d, psi, C, i, j) result (return)
    integer:: d
    integer,pointer,intent(in),contiguous:: i(:), j(:) 
    real(real_t):: return(size(i), size(j))
    real(real_t),pointer,intent(in),          &
      contiguous:: psi(:,:), C(:,:)           
    return = (                                &
      F(                                      &
        pi(d, psi, i,   j),                   &
        pi(d, psi, i+1, j),                   &
        pi(d,   C, i+h, j)                    &
      ) -                                     &
      F(                                      &
        pi(d, psi, i-1, j),                   &
        pi(d, psi, i,   j),                   &
        pi(d, C, i-h, j)                      &
      )                                       &
    )
  end function
!listing11
  subroutine donorcell_op_2D(psi, n, C, i, j)  
    class(arrvec_t),pointer:: psi, C
    integer,intent(in):: n
    integer,pointer,intent(in),contiguous:: i(:), j(:) 
    
    psi%at(n+1)%p%a(i,j) = psi%at(n)%p%a(i,j) &
      - donorcell(                            &
        0, psi%at(n)%p%a, C%at(0)%p%a, i,j    &
      )                                       &
      - donorcell(                            &
        1, psi%at(n)%p%a, C%at(1)%p%a, j,i    &
      )      
  end subroutine
!listing12
end module
!listing13
module donorcell_2D_m
  use donorcell_m
  use solver_2D_m
  implicit none
  
  type,extends(solver_2D_t):: donorcell_2D_t
    contains
    procedure:: ctor => donorcell_2D_ctor
    procedure:: dtor => donorcell_2D_dtor
    procedure:: advop => donorcell_2D_advop
  end type

  contains
  
  subroutine donorcell_2D_ctor(              &
    this, bcx, bcy, nx, ny                   &
  )
    class(donorcell_2D_t):: this
    class(bcd_t),intent(in),pointer:: bcx, bcy
    integer,intent(in):: nx, ny
    call solver_2D_ctor(this, bcx,bcy, nx,ny, 1)
  end subroutine

  subroutine donorcell_2D_advop(this)
    class(donorcell_2D_t):: this
    call donorcell_op_2D(                     &
      this%psi, this%n, this%C, this%i, this%j&
    )
  end subroutine

  subroutine donorcell_2D_dtor(this)
    class(donorcell_2D_t):: this
    call solver_2D_dtor(this)
  end subroutine
end module
!listing14
module mpdata_m
  use arrvec_m
  use arakawa_c_m
  use pi_m
  implicit none
  contains 
!listing15
  function frac(nom, den) result (return)
    real(real_t),intent(in):: nom(:,:), den(:,:)
    real(real_t):: return(size(nom,1), size(nom,2))
    where (den > 0)
      return = nom / den
    elsewhere
      return = 0
    end where
  end function
!listing16
  function A(d, psi, i, j) result (return)
    integer:: d
    real(real_t),pointer,intent(in),contiguous:: psi(:,:)
    integer,pointer,contiguous:: i(:), j(:) 
    real(real_t):: return(size(i), size(j))
    return = frac(                           &
      pi(d, psi, i+1, j) - pi(d, psi, i, j), &
      pi(d, psi, i+1, j) + pi(d, psi, i, j)  &
    )  
  end function
!listing17
  function B(d, psi, i, j) result (return)
    integer:: d
    real(real_t),pointer,intent(in),contiguous:: psi(:,:) 
    integer,pointer,contiguous:: i(:), j(:) 
    real(real_t):: return(size(i), size(j))
    return = frac(                           &
      pi(d, psi, i+1, j+1)                   &
    + pi(d, psi, i,   j+1)                   &
    - pi(d, psi, i+1, j-1)                   &
    - pi(d, psi, i,   j-1),                  &
      pi(d, psi, i+1, j+1)                   &
    + pi(d, psi, i,   j+1)                   &
    + pi(d, psi, i+1, j-1)                   &
    + pi(d, psi, i,   j-1)                   &
    ) / 2
  end function
!listing18
  function C_bar(d, C, i, j) result (return)
    integer:: d
    real(real_t),pointer,intent(in),contiguous:: C(:,:) 
    integer,pointer,contiguous:: i(:), j(:) 
    real(real_t):: return(size(i), size(j))

    return = (                               &
      pi(d, C, i+1, j+h) +                   &
      pi(d, C, i,   j+h) +                   &
      pi(d, C, i+1, j-h) +                   &
      pi(d, C, i,   j-h)                     &
    ) / 4               
  end function
!listing19
  function antidiff_2D(d, psi, i, j, C) result (return)
    integer:: d
    integer, pointer, contiguous :: i(:), j(:) 
    real(real_t):: return(size(i), size(j))
    real(real_t),pointer,intent(in),contiguous:: psi(:,:) 
    class(arrvec_t),pointer:: C
    return =                                  &
      abs(pi(d, C%at(d)%p%a, i+h, j))         &
      * (1 - abs(pi(d, C%at(d)%p%a, i+h, j))) &
      * A(d, psi, i, j)                       &
      - pi(d, C%at(d)%p%a, i+h, j)            &
      * C_bar(d, C%at(d-1)%p%a, i, j)         &
      * B(d, psi, i, j)
  end function
!listing20
end module
!listing21
module mpdata_2D_m
  use solver_2D_m
  use mpdata_m
  use donorcell_m
  use halo_m
  implicit none
  
  type, extends(solver_2D_t):: mpdata_2D_t
    integer :: n_iters, n_tmp
    integer, pointer, contiguous :: im(:), jm(:)
    class(arrvec_t), pointer:: tmp(:) 
    contains
    procedure:: ctor => mpdata_2D_ctor
    procedure:: dtor => mpdata_2D_dtor
    procedure:: advop => mpdata_2D_advop
  end type

  contains

  subroutine mpdata_2D_ctor(this, n_iters, bcx, bcy, nx, ny)
    class(mpdata_2D_t) :: this
    class(bcd_t), pointer :: bcx, bcy
    integer, intent(in) :: n_iters, nx, ny
    integer:: hlo, c

    hlo = 1
    this%n_iters = n_iters

    call solver_2D_ctor(this, bcx, bcy, nx, ny, hlo)

    this%n_tmp = 1
    if (n_iters > 2) this%n_tmp = 2
    allocate(this%tmp(0:this%n_tmp)) 

    do c=0, this%n_tmp - 1
      call this%tmp(c)%ctor(2)
      call this%tmp(c)%init(0, this%i // h, this%j // hlo)
      call this%tmp(c)%init(1, this%i // hlo, this%j // h)
    end do

    allocate(this%im(0 : nx))
    this%im = (/ (c, c=this%i(0) - 1, this%i(nx-1)) /)

    allocate(this%jm(0 : ny))
    this%jm = (/ (c, c=this%j(0) - 1, this%j(ny-1)) /)
  end subroutine

  subroutine mpdata_2D_dtor(this)
    class(mpdata_2D_t):: this
    integer:: c
    do c=0, this%n_tmp-1 
      call this%tmp(c)%dtor()
    end do
    call solver_2D_dtor(this)
    deallocate(this%im, this%jm)
  end subroutine

  subroutine mpdata_2D_advop(this)
    class(mpdata_2D_t):: this
    integer:: step

    do step=0, this%n_iters-1
      if (step == 0) then
        call donorcell_op_2D(this%psi, this%n, &
          this%C, this%i, this%j)
      else
        call this%cycle()
        call this%bcx%fill_halos(this%psi%at(this%n)%p%a, this%j // this%hlo)
        call this%bcy%fill_halos(this%psi%at(this%n)%p%a, this%i // this%hlo)

        block
          ! chosing input/output for antidiff. C
          class(arrvec_t), pointer :: C_corr, C_unco
          if (step == 1) then
            C_unco => this%C
            C_corr => this%tmp(0)
          else if (mod(step, 2) == 1) then
            C_unco => this%tmp(1) ! odd step
            C_corr => this%tmp(0) ! even step
          else
            C_unco => this%tmp(0) ! odd step
            C_corr => this%tmp(1) ! even step
          end if

          ! calculating the antidiffusive velo
          C_corr%at( 0 )%p%a( this%im+h, this%j ) &
            = antidiff_2D(0,                      &
              this%psi%at( this%n )%p%a,          & 
              this%im, this%j, C_unco             &
          )
          call this%bcy%fill_halos(C_corr%at(0)%p%a, this%i // h)

          C_corr%at( 1 )%p%a( this%i, this%jm+h ) &
            = antidiff_2D(1,                      &
              this%psi%at( this%n )%p%a,          &
              this%jm, this%i, C_unco             &
          )
          call this%bcx%fill_halos(C_corr%at(1)%p%a, this%j // h)

          ! donor-cell step
          call donorcell_op_2D(this%psi, this%n, C_corr, this%i, this%j)
        end block
      end if
    end do
  end subroutine
end module
!listing22
