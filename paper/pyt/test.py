import sys
from numpy import loadtxt
from numpy.testing import assert_almost_equal
from listings import *

if (len(sys.argv) != (9 + 1)) : 
  raise Exception('expecting 9 arguments (nx, ny, Cx, Cy, nt, it, f.in, f.out, dec)')

nx = int(sys.argv[1])
ny = int(sys.argv[2])
Cx = float(sys.argv[3])
Cy = float(sys.argv[4])
nt = int(sys.argv[5])
it = int(sys.argv[6])
fin = sys.argv[7]
fout = sys.argv[8]
dec = int(sys.argv[9])

slv = Solver_2D(Mpdata, it, Cyclic, Cyclic, nx, ny)
slv.state()[:] = loadtxt(fin)
slv.Cx()[:] = Cx
slv.Cy()[:] = Cy
slv.solve(nt)

assert_almost_equal(slv.state(), loadtxt(fout), dec)
