import numpy as np 
import subprocess
from scipy.io.netcdf import netcdf_file
import time

for lang in ('cpp','for','pyt') :
  # simulation parameters (via the netCDF global attributes)
  ncfile = 'data-' + lang + '.nc'
  nc = netcdf_file(ncfile, mode='w')
  nc.np = 1
  nc.nx = 10
  nc.ny = 10
  nc.no = 5
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
  
  t0 = time.time()
  cmd = ('./egu2012-' + lang,) if lang != 'pyt' else ('python', './egu2012.py')
  subprocess.check_call(cmd + (ncfile,))
  print time.time() - t0
  
  nc = netcdf_file(ncfile, mode='r')

  # plotting the result
  for t in range(nc.nt / nc.no + 1):
    print nc.variables['psi'][t,:,:]
