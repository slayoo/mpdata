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
    from numba import jit, float64, int_, void
    float64_2d=float64[:,:]
    print "Numba or Numbapro imported"
except ImportError:
    def void(*args):
        def dont_decorate(fn):
            return fn
        return dont_decorate
    float64 = void
    float64_2d = void
    int_ = void
    def jit(cl):
        return cl
    print "Not using Numba nor Numbapro"


import pdb


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
  @float64()
  def state(self):
    #pdb.set_trace()
    return self.psi[self.n][self.i]

  @void()
  def advop(self):
    #pdb.set_trace()                                                                       
    donorcell_op(
      self.psi, self.n,
      self.C, self.i
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
  

#listing09
def f(psi_l, psi_r, C):
  #pdb.set_trace()
  return (
    (C + abs(C)) * psi_l + 
    (C - abs(C)) * psi_r
  ) / 2
#listing10
def donorcell(psi, C, i):
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
def donorcell_op(psi, n, C, i):
  #pdb.set_trace()
  psi[n+1][i] = (psi[n][i] 
    - donorcell(psi[n], C, i)
     )

def main(psi_in, cx, nx=3, nt=1):
  slv = Donorcell(nx,1)
  slv.state()[:] = psi_in
  #pdb.set_trace()
  slv.C[:] = cx
  #pdb.set_trace()

  slv.solve(nt)
  print "psi_out", slv.state()[:]

Psi_in = numpy.array([0, 1, 0])
Cx = 1

main(Psi_in, Cx)
