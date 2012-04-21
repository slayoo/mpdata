import numpy as np 
import subprocess
from scipy.io.netcdf import netcdf_file
import time
from pylab import *

times = {}
sizes = (10,20,40,80,160,300)

for lang in ('cpp','for','pyt') :
  times[lang] = []
  for n in sizes :
    # simulation parameters (via the netCDF global attributes)
    ncfile = 'data-' + lang + '.nc'
    nc = netcdf_file(ncfile, mode='w')
    nc.np = 1
    nc.nx = n
    nc.ny = n
    nc.no = 1000
    nc.nt = 100
    nc.Cx = .5
    nc.Cy = .5

    nc.createDimension('t', 0) # unlimited dim
    nc.createDimension('x', nc.nx)
    nc.createDimension('y', nc.ny)
    psi = nc.createVariable('psi', 'float32', ('t','x','y'))
    psi[0,0,0] = 1

    # closing the file, running the solver and reopening the file
    nc.close()
  
    cmd = ('./egu2012-' + lang,) if lang != 'pyt' else ('python', './egu2012.py')
    t0 = time.time()
    subprocess.check_call(cmd + (ncfile,))
    times[lang].append(time.time() - t0)
  
    #nc = netcdf_file(ncfile, mode='r')
    # plotting the result
    #for t in range(nc.nt / nc.no + 1):
    #  print lang, nc.variables['psi'][t,:,:]

sizes =np.array([10,20,40,80,160,300])**2
figure(1)
plot(
  sizes, times['cpp'], 
  sizes, times['for'], 
  sizes, times['pyt']
)
legend(['C++','FORTRAN','Python'])
show()
