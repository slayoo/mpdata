module adv
  implicit none

  type adv_t
    integer :: ntlev = 2, halo = 1
    contains
    procedure :: left_halo
    procedure :: rght_halo
    procedure :: left_edge
    procedure :: rght_edge
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
  implicit none

  type, extends(adv_t) :: adv_upstream_t
  end type 
end module

program egu2012
  use netcdf
  use adv_upstream
  implicit none

  ! Fortran does not allow an array po pointers but it allows an array of derived types
  type :: arr_t
    real, dimension(:,:), pointer :: arr
  end type 

  character (len=256) :: fname
  integer :: ncid, stat, nt, nx, ny, np, no, i, t
  type(arr_t), allocatable :: psi(:)
  type(arr_t), pointer :: psi_ptr(:)
  
  class(adv_t), allocatable :: a
  real, dimension(:,:), pointer :: tmp

  ! handling command-line argument
  if (command_argument_count() /= 1) stop "ERROR - expecting one argument - a netCDF file name"
  call get_command_argument(1,fname)
  
  ! opening the netCDF file
  call check(nf90_open(fname, NF90_NOWRITE, ncid) )

  ! reading options
  call check(nf90_get_att(ncid, NF90_GLOBAL, "nx", nx))
  call check(nf90_get_att(ncid, NF90_GLOBAL, "ny", ny))
  call check(nf90_get_att(ncid, NF90_GLOBAL, "np", np))
  call check(nf90_get_att(ncid, NF90_GLOBAL, "no", no))
  call check(nf90_get_att(ncid, NF90_GLOBAL, "nt", nt))

  ! instantiating the advection operator
  allocate(adv_upstream_t::a)

  ! allocating memory
  allocate(psi(0:a%ntlev-1))
  allocate(psi_ptr(0:a%ntlev-1))
  do i = 0, a%ntlev - 1
    allocate(psi(i)%arr(-a%halo:nx-1+a%halo, -a%halo:ny-1+a%halo), stat=stat)
    if (stat /= 0) stop "ERROR - failed to allocate memory"
    psi_ptr(i)%arr => psi(i)%arr
  end do

  ! integration loop
  do t = 1, nt 
    ! filling halos (TODO: that's just one dimension as of now)
    tmp => a%left_halo(psi_ptr(0)%arr) 
    tmp = a%rght_edge(psi_ptr(0)%arr)
    tmp => a%left_edge(psi_ptr(0)%arr)
    tmp = a%rght_halo(psi_ptr(0)%arr)
     
    ! advecting in each dimension
    
    ! outputting
    if (modulo(t, no) == 0) print*, t, psi_ptr(1)%arr(0,0)
    
    ! cycling the pointers
    psi_ptr(0)%arr => psi(modulo(t, 2))%arr
    psi_ptr(1)%arr => psi(modulo(t+1, 2))%arr
    
  end do

  ! deallocating memory
  do i = 0, a%ntlev - 1
    deallocate(psi(i)%arr,stat=stat)
    if (stat /= 0) stop "ERROR - failed to deallocate memory"
  end do
  deallocate(psi,stat=stat)
  if (stat /= 0) stop "ERROR - failed to deallocate memory"

  ! closing the netCDF file
  call check(nf90_close(ncid))

contains
  subroutine check(status)
    integer, intent(in) :: status
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "ERROR - while accessing the netCDF file"
    end if
  end subroutine check 

end program egu2012
