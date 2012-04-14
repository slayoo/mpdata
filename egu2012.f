program egu2012
  use netcdf
  implicit none

  ! Fortran does not allow an array of pointers but it allows an array of derived types
  type :: arr_t
    real, dimension(:,:), pointer :: arr
  end type 

  character (len=256) :: fname
  integer :: ncid, stat, nt, nx, ny, np, no, i, t
  type(arr_t), allocatable :: psi(:)
  type(arr_t), pointer :: psi_ptr(:)
  
  integer, parameter :: ntlev = 2, halo = 1! TODO: adv::ntlev, ...

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

  ! allocating memory
  allocate(psi(0:ntlev-1))
  allocate(psi_ptr(0:ntlev-1))
  do i = 0, ntlev - 1
    allocate(psi(i)%arr(-halo:nx-1+halo, -halo:ny-1+halo),stat=stat)
    if (stat /= 0) stop "ERROR - failed to allocate memory"
    psi_ptr(i)%arr => psi(i)%arr
  end do

  ! integration loop
  do t = 1, nt 
    ! filling haloes
     
    ! advecting in each dimension
    
    ! outputting
    if (modulo(t, no) == 0) print*, t, psi(1)%arr(0,0)
    
    ! cycling the pointers
    psi_ptr(0)%arr => psi(modulo(t, 2))%arr
    psi_ptr(1)%arr => psi(modulo(t, 2))%arr
    
  end do

  ! deallocating memory
  do i = 0, ntlev - 1
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
