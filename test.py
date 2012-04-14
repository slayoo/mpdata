import numpy as np 
import subprocess
from scipy.io.netcdf import netcdf_file

# simulation parameters (via the netCDF global attributes)
ncfile = 'data.nc'
nc = netcdf_file(ncfile, mode='w')
nc.np = 1
nc.nx = 3
nc.ny = 4
nc.no = 1
nc.nt = 10
nc.Cx = .5
nc.Cy = .5

nc.createDimension('t', 0) # unlimited dim
nc.createDimension('x', nc.nx)
nc.createDimension('y', nc.ny)
psi = nc.createVariable('psi', 'float32', ('t','x','y'))
psi[0,0,0] = 1

# closing the file, running the solver and reopening the file
nc.close()
subprocess.check_call(('./egu2012', ncfile))
nc = netcdf_file(ncfile, mode='r')

# plotting the result
for t in range(nc.nt / nc.no + 1):
  print nc.variables['psi'][t,:,:]
  # TODO!
