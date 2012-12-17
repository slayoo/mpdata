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
except ImportError:
  pass
import numpy
import pdb
#listing03
class Shift():
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
one = Shift(1,1) 
hlf = Shift(0,1)
#listing05
def pi(d, *idx): 
  return (idx[d], idx[d-1])
#listing06
class Solver_2D(object):
  # ctor-like method
  def __init__(self, bcx, bcy, nx, ny, hlo):
    self.n = 0
    self.i = slice(hlo, nx + hlo)
    self.j = slice(hlo, ny + hlo)

    self.bcx = bcx(0, self.i, hlo)
    self.bcy = bcy(1, self.j, hlo)

    self.psi = (
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t),
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t) 
    )
    self.C = (
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t),
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t)
    )

  # accessor methods
  def state(self):
    return self.psi[self.n][self.i, self.j]

  # helper methods invoked by solve()
  def courant(self,d):
    return self.C[d][self.i, self.j]

  def cycle(self):
    self.n  = (self.n + 1) % 2 - 2

  def xchng(self, arr):
    self.bcx.fill_halos(arr)
    self.bcy.fill_halos(arr)

   # integration logic
  def solve(self, nt):
    for t in range(nt):
      #pdb.set_trace()
      self.xchng(self.C[0])
      self.xchng(self.C[1])
      self.xchng(self.psi[self.n])
      self.advop() 
      self.cycle()
  
#listing07
class Cyclic(object):
  # ctor
  def __init__(self, d, i, hlo): 
    self.left_halo = pi(d, 
      slice(i.start-hlo, i.start), slice(None)
    )
    self.rght_edge = pi(d, 
      slice(i.stop-hlo, i.stop),   slice(None)
    )
    self.rght_halo = pi(d, 
      slice(i.stop, i.stop+hlo),   slice(None)
    )
    self.left_edge = pi(d, 
      slice(i.start, i.start+hlo), slice(None)
    )

  # method invoked by the solver
  def fill_halos(self, psi):
    psi[self.left_halo] = psi[self.rght_edge]
    psi[self.rght_halo] = psi[self.left_edge]

#listing08
def f(psi_l, psi_r, C):
  return (
    (C + abs(C)) * psi_l + 
    (C - abs(C)) * psi_r
  ) / 2
#listing09
def donorcell(d, psi, C, i, j):
  #pdb.set_trace()
  return (
    f(
      psi[pi(d, i,     j)], 
      psi[pi(d, i+one, j)], 
        C[pi(d, i+hlf, j)]
    ) - 
    f(
      psi[pi(d, i-one, j)], 
      psi[pi(d, i,     j)], 
        C[pi(d, i-hlf, j)]
    ) 
  )
#listing10
def donorcell_op_2D(psi, n, C, i, j):
  psi[n+1][i,j] = (psi[n][i,j] 
    - donorcell(0, psi[n], C[0], i, j)
    - donorcell(1, psi[n], C[1], j, i)
  )

class Donorcell_2D(Solver_2D):
  def __init__(self, bcx, bcy, nx, ny):
    Solver_2D.__init__(self, bcx, bcy, nx, ny, 1)

  def advop(self):
    donorcell_op_2D(self.psi, self.n, self.C, self.i, self.j)

#listing11
def frac(nom, den):
  return numpy.where(den > 0, nom/den, 0)

#listing12
def a_op(d, psi, i, j):
  return frac(
    psi[pi(d, i+one, j)] - psi[pi(d, i, j)],
    psi[pi(d, i+one, j)] + psi[pi(d, i, j)]
  )
#listing13
def b_op(d, psi, i, j):
  return frac( 
      psi[pi(d, i+one, j+one)] 
    + psi[pi(d, i,     j+one)] 
    - psi[pi(d, i+one, j-one)] 
    - psi[pi(d, i,     j-one)],
      psi[pi(d, i+one, j+one)] 
    + psi[pi(d, i,     j+one)] 
    + psi[pi(d, i+one, j-one)] 
    + psi[pi(d, i,     j-one)]
  ) / 2
#listing14
def C_bar(d, C, i, j):
  return (
    C[pi(d, i+one, j+hlf)] + 
    C[pi(d, i,     j+hlf)] +
    C[pi(d, i+one, j-hlf)] + 
    C[pi(d, i,     j-hlf)] 
  ) / 4
#listing14
def antidiff_2D(d, psi, i, j, C):
  return (
    abs(C[d][pi(d, i+hlf, j)]) 
    * (1 - abs(C[d][pi(d, i+hlf, j)])) 
    * a_op(d, psi, i, j)
    - C[d][pi(d, i+hlf, j)] 
    * C_bar(d, C[d-1], i, j)
    * b_op(d, psi, i, j)
  )
#listing15
class Mpdata_2D(Solver_2D):
  def __init__(self, n_iters, bcx, bcy, nx, ny):
    hlo = 1
    Solver_2D.__init__(self, bcx, bcy, nx, ny, hlo)

    self.n_iters = n_iters
  
    self.tmp = [(
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t),
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t)
    )]
    
    if n_iters > 2:
      self.tmp.append((
        numpy.empty(( nx+2*hlo, ny+2*hlo), real_t),
        numpy.empty(( nx+2*hlo, ny+2*hlo), real_t)
      ))

  def advop(self):
    for step in range(self.n_iters):
      if step == 0:
        donorcell_op_2D(self.psi, self.n, self.C, self.i, self.j)
      else:
        self.cycle()
        self.xchng(self.psi[self.n])
        if step == 1:
          C_unco, C_corr = self.C, self.tmp[0]
        elif step % 2:
          C_unco, C_corr = self.tmp[1], self.tmp[0]
        else:
          C_unco, C_corr = self.tmp[0], self.tmp[1]

        im = self.i
        jm = self.j

        C_corr[0][im+hlf, self.j] = (
          antidiff_2D(0, self.psi[self.n], im, self.j, C_unco)) 
        self.xchng(C_corr[0])

        C_corr[1][self.i, jm+hlf] = (
          antidiff_2D(1, self.psi[self.n], jm, self.i, C_unco)) 
        self.xchng(C_corr[1])

        pdb.set_trace()
        donorcell_op_2D(self.psi, self.n, C_corr, self.i, self.j)
