#listing00
# code licensed under the terms of GNU GPL v3
# copyright holder: University of Warsaw
#listing01
real_t = 'float64'
#listing02
#try:
#  import numpypy
#  from _numpypy.pypy import set_invalidation
#  set_invalidation(False)
#except ImportError:
#  pass
cimport numpy
import numpy
try: 
  numpy.seterr(all='ignore')
except AttributeError:
  pass
#listing03
class shift():
  def __init__(self, plus, mnus):
    self.plus = plus
    self.mnus = mnus
  def __radd__(self, arg): 
    return type(arg)(
      arg.start + self.plus, 
      arg.stop  + self.plus
    )
  def __rsub__(self, arg): 
    return type(arg)(
      arg.start - self.mnus, 
      arg.stop  - self.mnus
    )
#listing04
one = shift(1,1) 
hlf = shift(0,1)
#listing05
def ext(r, n):
  if (type(n) == int) & (n == 1): 
    n = one
  return slice(
    (r - n).start, 
    (r + n).stop
  )
#listing06
def pi(d, *idx): 
  return (idx[d], idx[d-1])
#listing07
class solver(object):
  # ctor-like method
  def __init__(self, bcx, bcy, nx, ny, hlo):
    self.n = 0
    self.hlo = hlo
    self.i = slice(hlo, nx + hlo)
    self.j = slice(hlo, ny + hlo)

    self.bcx = bcx(0, self.i, hlo)
    self.bcy = bcy(1, self.j, hlo)

    self.psi = (
      numpy.empty((
        ext(self.i, self.hlo).stop, 
        ext(self.j, self.hlo).stop
      ), real_t),
      numpy.empty((
        ext(self.i, self.hlo).stop, 
        ext(self.j, self.hlo).stop
      ), real_t)
    )

    self.C = (
      numpy.empty((
        ext(self.i, hlf).stop, 
        ext(self.j, self.hlo).stop
      ), real_t),
      numpy.empty((
        ext(self.i, self.hlo).stop, 
        ext(self.j, hlf).stop
      ), real_t)
    )

  # accessor methods
  def state(self):
    return self.psi[self.n][self.i, self.j]

  # helper methods invoked by solve()
  def courant(self,d):
    return self.C[d][:]

  def cycle(self):
    self.n  = (self.n + 1) % 2 - 2

   # integration logic
  def solve(self, nt):
    for t in range(nt):
      #print "t = ", t
      self.bcx.fill_halos(
        self.psi[self.n], ext(self.j, self.hlo)
      )
      self.bcy.fill_halos(
        self.psi[self.n], ext(self.i, self.hlo)
      )
      self.advop() 
      self.cycle()
  
#listing08
class cyclic(object):
  # ctor
  def __init__(self, d, i, hlo): 
    self.d = d
    self.left_halo = slice(i.start-hlo, i.start    )
    self.rght_edge = slice(i.stop -hlo, i.stop     )
    self.rght_halo = slice(i.stop,      i.stop +hlo)
    self.left_edge = slice(i.start,     i.start+hlo)

  # method invoked by the solver
  def fill_halos(self, psi, j):
    psi[pi(self.d, self.left_halo, j)] = (
      psi[pi(self.d, self.rght_edge, j)]
    )
    psi[pi(self.d, self.rght_halo, j)] = (
      psi[pi(self.d, self.left_edge, j)]
    )

#listing09
def f(double psi_l, double psi_r, double C):
  return (
    (C + abs(C)) * psi_l + 
    (C - abs(C)) * psi_r
  ) / 2

#listing10
def donorcell_x(numpy.ndarray[double, ndim=2] psi, numpy.ndarray[double, ndim=2] C, 
    				      int i, int j):
  return (
    f(psi[i,   j],
      psi[i+1, j], 
      C[i, j]
      ) - 
    f(
      psi[i-1, j], 
      psi[i,   j], 
        C[i-1, j]
      ) 
    )

def donorcell_y(numpy.ndarray[double, ndim=2] psi, numpy.ndarray[double, ndim=2] C, 
    int i, int j):
  return (
    f(
      psi[i, j], 
      psi[i, j+1], 
        C[i, j]
    ) - 
    f(
      psi[i, j-1], 
      psi[i, j], 
        C[i, j-1]
    ) 
  )

#listing11
def donorcell_op(tuple psi, int n, tuple C, slice i_sl, slice j_sl):
  i_range = range(i_sl.start, i_sl.stop)
  j_range = range(j_sl.start, j_sl.stop)
  cdef numpy.ndarray[double, ndim=2] psi_f = psi[n+1]
  cdef numpy.ndarray[double, ndim=2] psi_i = psi[n]
  cdef numpy.ndarray[double, ndim=2] C_0 = C[0]
  cdef numpy.ndarray[double, ndim=2] C_1 = C[1]
  cdef unsigned int i, j

  for i in i_range:
    for j in j_range:
      psi_f[i,j] = (psi_i[i,j] 
    - donorcell_x(psi_i, C_0, i, j)
    - donorcell_y(psi_i, C_1, i, j)
  )

#listing12
class solver_donorcell(solver):
  def __init__(self, bcx, bcy, nx, ny):
    solver.__init__(self, bcx, bcy, nx, ny, 1)

  def advop(self):
    donorcell_op(
      self.psi, self.n, 
      self.C, self.i, self.j
    )
