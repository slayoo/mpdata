from listings import *
import numpy 

nx, ny = 3, 3
Cx, Cy = 0, 1
psi_in = numpy.array([[0,0,0],[0,1,0],[0,0,0]])

slv = Mpdata_2D(2, Cyclic, Cyclic, nx, ny)
slv.state()[:] = psi_in
slv.courant(0)[:] = Cx
slv.courant(1)[:] = Cy
slv.solve(1)

print slv.state()

#if (abs(slv.state() - read_file(fout, nx, ny)) >= .5 * pow(10, -dec)).any(): 
#  print slv.state()
#  raise Exception()
