from listings import *
import sys

def read_file(fname, nx, ny):
  tmp = numpy.empty((nx, ny), real_t)
  with open(fname, 'r') as f:
    for x in range(0, nx):
      tmp[x,:] = numpy.fromstring(f.readline(), dtype=real_t, sep='\t')
  return tmp

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
slv.state()[:] = read_file(fin, nx, ny)
slv.Cx()[:] = Cx
slv.Cy()[:] = Cy
slv.solve(nt)

if (abs(slv.state() - read_file(fout, nx, ny)) >= .5 * pow(10, -dec)).any(): 
  raise Exception()
