module hack
  ! Fortran does not allow an array of pointers but it allows an array of derived types
  type :: arr_t
    real, dimension(:,:), pointer :: X
  end type 
end module

module adv
  use hack
  implicit none

  type, abstract :: adv_t
    integer, public :: ntlev = 2, halo = 1, mh = 0, ph = 1
    contains   
    procedure, public :: rng_psi, rng_vel 
    procedure, public :: left_halo, rght_halo, left_edge, rght_edge
    procedure(op_proto), public, deferred :: op
  end type 

  abstract interface 
    subroutine op_proto(this, psi, vel, i, j, n, d)
      import :: adv_t
      import :: arr_t
      class(adv_t), intent(in) :: this
      type(arr_t), intent(inout), pointer :: psi(:)
      type(arr_t), intent(in), pointer :: vel(:)
      integer, intent(in), dimension(:) :: i, j
      integer, intent(in) :: n, d
    end subroutine
  end interface 

  contains

  function left_halo(this, arr, d)
    class(adv_t), intent(in) :: this
    real, intent(in), dimension(:,:), pointer :: arr
    integer, intent(in) :: d
    real, dimension(:,:), pointer :: left_halo
    select case (d)
      case (1)
        left_halo => arr(-this%halo : -1, :)
      case (2)
        left_halo => arr(:, -this%halo : -1)
      case default 
        stop
    end select
  end function 

  function rght_halo(this, arr, d)
    class(adv_t), intent(in) :: this
    real, intent(in), dimension(:,:), pointer :: arr
    integer, intent(in) :: d
    real, dimension(:,:), pointer :: rght_halo
    select case (d)
      case (1)
        rght_halo => arr(size(arr,d)-2*this%halo : size(arr,d)-this%halo-1, :)
      case (2)
        rght_halo => arr(:, size(arr,d)-2*this%halo : size(arr,d)-this%halo-1)
      case default 
        stop
    end select
  end function 

  function left_edge(this, arr, d)
    class(adv_t), intent(in) :: this
    real, intent(in), dimension(:,:), pointer :: arr
    integer, intent(in) :: d
    real, dimension(:,:), pointer :: left_edge
    select case (d)
      case (1)
        left_edge => arr(0 : this%halo - 1, :)
      case (2)
        left_edge => arr(:, 0 : this%halo - 1)
      case default 
        stop
    end select
  end function 

  function rght_edge(this, arr, d)
    class(adv_t), intent(in) :: this
    real, intent(in), dimension(:,:), pointer :: arr
    integer, intent(in) :: d
    real, dimension(:,:), pointer :: rght_edge
    select case (d)
      case (1)
        rght_edge => arr(size(arr,d)-3*this%halo : size(arr,d)-2*this%halo-1, :)
      case (2)
        rght_edge => arr(:, size(arr,d)-3*this%halo : size(arr,d)-2*this%halo-1)
      case default 
        stop
    end select
  end function 

end module

module adv_upstream
  use adv
  use hack
  implicit none

  type, extends(adv_t) :: adv_upstream_t
    contains
    procedure, public :: op
    !procedure, private :: F
  end type 

  contains 

  ! eq. (3 a-d) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480)
  ! statement functions are obsolete and do not support array-valued arguments
  ! vector-indexing results in copying the argument, hence restoring to preprocessor :(
#define F(psi_l, psi_r, U) (max(0.,U) * psi_l + min(0.,U) * psi_r)

  ! eq. (2) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) 
  subroutine op(this, psi, vel, i, j, n, d)
    class(adv_upstream_t), intent(in) :: this
    type(arr_t), intent(inout), pointer :: psi(:)
    type(arr_t), intent(in), pointer :: vel(:)
    integer, intent(in), dimension(:) :: i, j
    integer, intent(in) :: n, d
    
    select case (d)
      case (1)
        psi( n+1 )%X( i,j ) = psi( n+1 )%X( i,j ) - (                            & 
          F(psi( n )%X( i,  j ), psi( n )%X( i+1,j ), vel(0)%X( i,j+this%ph )) - &
          F(psi( n )%X( i-1,j ), psi( n )%X( i,  j ), vel(0)%X( i,j-this%mh ))   &
        )
      case (2)
        psi( n+1 )%X( i,j ) = psi( n+1 )%X( i,j ) - (                            & 
          F(psi( n )%X( i,j   ), psi( n )%X( i,j+1 ), vel(0)%X( i,j+this%ph )) - &
          F(psi( n )%X( i,j-1 ), psi( n )%X( i,j   ), vel(0)%X( i,j-this%mh ))   &
        )
      case default 
        stop
    end select
  end subroutine
end module

program egu2012
  use hack
  use netcdf
  use adv_upstream
  implicit none

  character (len=256) :: fname
  integer :: ncid, varid, stat, nt, nx, ny, np, no, i, j, p, t, d
  real :: Cx, Cy
  integer, parameter :: n = 0
  integer, dimension(:,:), allocatable, target :: ii, jj
  type(arr_t), allocatable :: psi(:)
  type(arr_t), pointer :: psi_ptr(:), vel(:)
  
  class(adv_upstream_t), allocatable :: a
  real, dimension(:,:), pointer :: tmp

  ! handling command-line argument
  if (command_argument_count() /= 1) stop "ERROR - expecting one argument - a netCDF file name"
  call get_command_argument(1,fname)
  
  ! opening the netCDF file
  stat = nf90_open(fname, NF90_WRITE, ncid) 

  ! reading options
  stat = nf90_get_att(ncid, NF90_GLOBAL, "nx", nx)
  stat = nf90_get_att(ncid, NF90_GLOBAL, "ny", ny)
  stat = nf90_get_att(ncid, NF90_GLOBAL, "np", np)
  stat = nf90_get_att(ncid, NF90_GLOBAL, "no", no)
  stat = nf90_get_att(ncid, NF90_GLOBAL, "nt", nt)
  stat = nf90_get_att(ncid, NF90_GLOBAL, "Cx", Cx)
  stat = nf90_get_att(ncid, NF90_GLOBAL, "Cy", Cy)

  ! instantiating the advection operator
  allocate(adv_upstream_t::a) 

  ! allocating memory for psi
  allocate(psi(0:a%ntlev-1), psi_ptr(0:a%ntlev-1)) 
  do i = 0, a%ntlev - 1
    allocate(psi(i)%X(-a%halo:nx-1+a%halo, -a%halo:ny-1+a%halo), stat=stat) !TODO: rng_psi
    psi_ptr(i)%X => psi(i)%X
  end do

  ! helper vars for vector indexing the arrays
  allocate(ii(np, 0:nx-1), jj(np, 0:ny-1))
  do p = 0, np - 1
    do i = p*nx/np, (p+1)*nx/np - 1 
      ii(p, i) = i
    end do
    do j = p*ny/np, (p+1)*ny/np - 1 
      jj(p, j) = j
    end do
  end do

  ! filling psi with data from netCDF
  stat = nf90_inq_varid(ncid, "psi", varid)
  stat = nf90_get_var(ncid, varid, psi(n)%X(0:nx-1, 0:ny-1))

  ! allocating memory for vel
  allocate(vel(0:1))
  allocate( &
    vel(0)%X(-a%halo:nx-1+a%halo, -a%halo:ny-1+a%halo),& !TODO: rng_vel 
    vel(1)%X(-a%halo:nx-1+a%halo, -a%halo:ny-1+a%halo) &
  )

  ! filling vel with data from netCDF
  vel(0)%X = Cx
  vel(1)%X = Cy

  ! integration loop
  do t = 1, nt 
    ! without the temporary variable GCC segfaults on compilation! (bug no. 52994)
    do d = 1, 2
      tmp => a%left_halo(psi_ptr(0)%X, d) 
      tmp = a%rght_edge(psi_ptr(0)%X, d)
      tmp => a%left_halo(psi_ptr(0)%X, d)
      tmp = a%rght_edge(psi_ptr(0)%X, d)
    end do
     
    ! advecting in each dimension
    psi_ptr(n+1)%X = psi_ptr(n)%X
!$OMP PARALLEL DO
    do p = 0, np - 1
      do d = 1, 2
        call a%op(psi_ptr, vel, ii(p,:), jj(p,:), n, d)
      end do
    end do
!$OMP END PARALLEL DO
    
    ! outputting
    if (modulo(t, no) == 0) then
      stat = nf90_put_var(ncid, varid, psi(n+1)%X(0:nx-1, 0:ny-1), start = (/1,1,1+t/no/))
    end if
    
    ! cycling the pointers
    psi_ptr(n)%X => psi(modulo(t, 2))%X
    psi_ptr(n+1)%X => psi(modulo(t+1, 2))%X
  end do

  ! cleanup: deallocating memory and closing the netcdf file
  do i = 0, a%ntlev - 1
    deallocate(psi(i)%X,stat=stat)
  end do
  deallocate(psi,stat=stat)
  stat = nf90_close(ncid)
end program egu2012
