module hack
  ! Fortran does not allow an array of pointers but it allows an array of derived types
  type :: arr_t
    real, dimension(:,:), pointer :: X
  end type 
end module

module adv
  implicit none

  type adv_t
    integer :: ntlev = 2, halo = 1
    contains
    procedure, public :: left_halo
    procedure, public :: rght_halo
    procedure, public :: left_edge
    procedure, public :: rght_edge
  end type 

  contains

  function left_halo(this, arr)
    class(adv_t), intent(in) :: this
    real, intent(in), dimension(:,:), pointer :: arr
    real, dimension(:,:), pointer :: left_halo
    left_halo => arr(-this%halo : -1, :)
  end function 

  function rght_halo(this, arr)
    class(adv_t), intent(in) :: this
    real, intent(in), dimension(:,:), pointer :: arr
    real, dimension(:,:), pointer :: rght_halo
    rght_halo => arr(size(arr,1)-2*this%halo : size(arr,1)-this%halo-1, :)
  end function 

  function left_edge(this, arr)
    class(adv_t), intent(in) :: this
    real, intent(in), dimension(:,:), pointer :: arr
    real, dimension(:,:), pointer :: left_edge
    left_edge => arr(0 : this%halo - 1, :)
  end function 

  function rght_edge(this, arr)
    class(adv_t), intent(in) :: this
    real, intent(in), dimension(:,:), pointer :: arr
    real, dimension(:,:), pointer :: rght_edge
    rght_edge => arr(size(arr,1)-3*this%halo : size(arr,1)-2*this%halo-1, :)
  end function 

end module

module adv_upstream
  use adv
  use hack
  implicit none

  type, extends(adv_t) :: adv_upstream_t
    contains
    procedure, public :: op
    procedure, private :: F
  end type 

  contains 

  ! eq. (3 a-d) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480)
  function F(this, psi_l, psi_r, U)
    class(adv_upstream_t), intent(in) :: this
    real, intent(in), dimension(:,:) :: psi_l
    real, intent(in), dimension(:,:) :: psi_r
    real, intent(in), dimension(:,:) :: u
    real, dimension(size(psi_l, 1), size(psi_l, 2)) :: F

    F = max(0.,U) * psi_l + min(0.,U) * psi_r 
  end function

  ! eq. (2) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) 
  subroutine op(this, psi, vel, i, j, n)
    class(adv_upstream_t), intent(in) :: this
    type(arr_t), intent(inout), pointer :: psi(:)
    type(arr_t), intent(in), pointer :: vel(:)
    integer, intent(in), dimension(:) :: i, j
    integer, intent(in) :: n

    psi( n+1 )%X( i,j ) = psi( n+1 )%X( i,j ) - (                           & 
      this%F(psi( n )%X( i,  j ), psi( n )%X( i+1,j ), vel(0)%X( i+1,j )) - & !TODO: phalf, dimension
      this%F(psi( n )%X( i-1,j ), psi( n )%X( i,  j ), vel(0)%X( i-0,j ))   & !TODO: mhalf, dimension
    )
  end subroutine
end module

program egu2012
  use hack
  use netcdf
  use adv_upstream
  implicit none

  character (len=256) :: fname
  integer :: ncid, stat, nt, nx, ny, np, no, i, t
  real :: Cx, Cy
  integer, parameter :: n = 0
  integer, dimension(:), allocatable :: ii, jj
  type(arr_t), allocatable :: psi(:)
  type(arr_t), pointer :: psi_ptr(:), vel(:)
  
  class(adv_upstream_t), allocatable :: a
  real, dimension(:,:), pointer :: tmp

  ! handling command-line argument
  if (command_argument_count() /= 1) stop "ERROR - expecting one argument - a netCDF file name"
  call get_command_argument(1,fname)
  
  ! opening the netCDF file
  stat = nf90_open(fname, NF90_NOWRITE, ncid) 

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

  ! allocating memory: psi
  allocate(psi(0:a%ntlev-1), psi_ptr(0:a%ntlev-1)) 
  do i = 0, a%ntlev - 1
    allocate(psi(i)%X(-a%halo:nx-1+a%halo, -a%halo:ny-1+a%halo), stat=stat) !TODO: rng_psi
    psi_ptr(i)%X => psi(i)%X
  end do

  ! allocating memory: vel
  allocate(vel(0:1))
  allocate( &
    vel(0)%X(-a%halo:nx-1+a%halo, -a%halo:ny-1+a%halo),& !TODO: rng_vel 
    vel(1)%X(-a%halo:nx-1+a%halo, -a%halo:ny-1+a%halo) &
  )
  vel(0)%X = 1!Cx
  vel(1)%X = Cy

  ! helpe... TODO
  allocate(ii(0:nx-1), jj(0:ny-1)) !TODO: zmienne dla OpenMP
  do i = 0, nx - 1 
    ii(i) = i
  end do
  do i = 0, ny - 1 
    jj(i) = i
  end do

  psi(0)%X(0,0) = 1 ! TODO: read from netCDF
  !print*, 0, psi_ptr(0)%X(ii,jj) ! TEMP

  ! integration loop
  do t = 1, nt 
    ! filling halos (TODO: that's just one dimension as of now)
    ! without the temporary variable GCC segfaults on compilation! (bug no. 52994)
    tmp => a%left_halo(psi_ptr(0)%X) 
    tmp = a%rght_edge(psi_ptr(0)%X)
    tmp => a%left_halo(psi_ptr(0)%X)
    tmp = a%rght_edge(psi_ptr(0)%X)
     
    ! advecting in each dimension
    
    psi_ptr(n+1)%X = psi_ptr(n)%X
    call a%op(psi_ptr, vel, ii, jj, n) ! TODO: second dimension
    call a%op(psi_ptr, vel, ii, jj, n) ! TODO: second dimension
    
    ! outputting
    !if (modulo(t, no) == 0) print*, t, psi_ptr(1)%X(ii,jj)
    
    ! cycling the pointers
    psi_ptr(n)%X => psi(modulo(t, 2))%X
    psi_ptr(n+1)%X => psi(modulo(t+1, 2))%X
    
  end do

  ! deallocating memory
  do i = 0, a%ntlev - 1
    deallocate(psi(i)%X,stat=stat)
  end do
  deallocate(psi,stat=stat)

  ! closing the netCDF file
  stat = nf90_close(ncid)

end program egu2012
