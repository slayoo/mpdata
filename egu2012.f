
program egu2012
  !use netcdf
  implicit none

  if (command_argument_count() .ne. 1) then
    print *, 'expecting one argument - a netCDF file name'
    STOP !call EXIT(1)
  end if

end program egu2012
