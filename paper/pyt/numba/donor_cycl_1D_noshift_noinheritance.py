#listing00
# code licensed under the terms of GNU GPL v3
# copyright holder: University of Warsaw
#listing01
real_t = 'float64'
#listing02
try:
  import numpypy
  from _numpypy.pypy import set_invalidation
  set_invalidation(False)
  print "Numpypy impported"
except ImportError:
  print "Not using numpypy"
  pass

import numpy
try: 
  numpy.seterr(all='ignore')
except AttributeError:
  pass

try:
    from numba import jit, float64, int_, void, autojit
    float64_1d=float64[:]
    float64_2d=float64[:,:]
    print "Numba or Numbapro imported"
except ImportError:
    def void(*args):
        def dont_decorate(fn):
            return fn
        return dont_decorate
    float64 = void
    float64_1d = void
    float64_2d = void
    int_ = void
    def jit(cl):
        return cl
    autojit = jit
    print "Not using Numba nor Numbapro"


import pdb

#listing09
#@autojit - done later explicitly
def f(psi_l, psi_r, C):
  #pdb.set_trace()
  return (
    (C + abs(C)) * psi_l + 
    (C - abs(C)) * psi_r
  ) / 2
#listing10
#@autojit - done later explicitly
def donorcell(psi, C, start, stop):
  i = slice(start,stop)
  i_pl_one = slice(i.start + 1, i.stop + 1)
  i_mn_one = slice(i.start - 1, i.stop - 1)
  i_pl_hlf = i
  i_mn_hlf = slice(i.start - 1, i.stop - 1)  
  #pdb.set_trace()
  return (
    f(psi[i],      psi[i_pl_one], C[i_pl_hlf]) - 
    f(psi[i_mn_one,], psi[i],     C[i_mn_hlf]) 
  )
#listing11
#@autojit - done later explicitly
def donorcell_op(psi, n, C, start, stop):
  i = slice(start,stop)
  #pdb.set_trace()
  psi[n+1][i] = (psi[n][i] 
    - donorcell(psi[n], C, start, stop)
     )

try:
  fast_f = jit(float64_1d, (float64_1d, float64_1d, float64_1d))(f)
  f = fast_f
  fast_donorcell = jit(float64_1d, (float64_1d, float64_1d, int_, int_))(donorcell)
  donorcell = fast_donorcell
  fast_donorcell_op = jit(void, (float64_2d, int_, float64_1d, int_, int_))(donorcell_op)
  donorcell_op = fast_donorcell_op
  print "Numba or Numbapro used"
except TypeError:
  print "Numba and Numbapro not working"
  
  

#listing07
@jit
class Donorcell(object):
  # ctor-like method
  def __init__(self, nx, hlo):
    self.n = 0
    self.hlo = hlo
    self.i = slice(hlo, nx + hlo)
   

    self.psi = (
      numpy.empty((
          self.i.stop - self.i.start + 2 * self.hlo), real_t),
      numpy.empty((
          self.i.stop - self.i.start + 2 * self.hlo), real_t)
      )

    self.C = numpy.empty((
        self.i.stop - self.i.start + self.hlo), real_t)
    #pdb.set_trace()

  # accessor methods
  @float64_1d()
  def state(self):
    #pdb.set_trace()
    return self.psi[self.n][self.i]

  @void()
  def advop(self):
    #pdb.set_trace()                                                                       
    donorcell_op(
      self.psi, self.n,
      self.C, self.i.start, self.i.stop
    )

  @void()
  def cycle(self):
    #pdb.set_trace()
    self.n  = (self.n + 1) % 2 - 2

  @void()
  def fill_halos(self):
    left_halo = slice(self.i.start-self.hlo, self.i.start    )
    rght_edge = slice(self.i.stop -self.hlo, self.i.stop     )
    rght_halo = slice(self.i.stop,           self.i.stop +self.hlo)
    left_edge = slice(self.i.start,          self.i.start+self.hlo)

    self.psi[self.n][left_halo] = self.psi[self.n][rght_edge]
    self.psi[self.n][rght_halo] = self.psi[self.n][left_edge]


   # integration logic
  @void(int_)
  def solve(self, nt):
    for t in range(nt):
      self.fill_halos()
      self.advop() 
      self.cycle()
  

def main(cx, nx=3, nt=1):
  psi_in = numpy.zeros(nx)
  psi_in[int(nx/4)] = 1
  slv = Donorcell(nx,1)
  slv.state()[:] = psi_in
  #pdb.set_trace()
  slv.C[:] = cx
  #pdb.set_trace()

  import time
  s = time.time()
  slv.solve(nt)
  e = time.time()
  print "psi_out", slv.state()[:].max()
  print "solved in", e - s


import sys
if len(sys.argv) == 2:
  main(1, nx=10000, nt=100000)
else:
  main(1)

