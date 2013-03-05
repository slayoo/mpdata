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

#listing03
@autojit
def sl_add(sl, int_num):
  return slice(sl.start + int_num, sl.stop + int_num)

#listing05
@autojit
def ext(r, lft, rgh):
  return slice(
    r.start - lft, 
    r.stop + rgh
  )
#listing06
@autojit
def pi(d, idx_0, idx_1):
  idx = [idx_0, idx_1]
  return (idx[d], idx[d-1]) 

#listing07
class Solver(object):
  # ctor-like method
  def __init__(self, bcx, bcy, nx, ny, hlo):
    self.n = 0
    self.hlo = hlo
    self.i = slice(hlo, nx + hlo)
    self.j = slice(hlo, ny + hlo)

    self.bcx = bcx(0, self.i.start, self.i.stop, hlo)
    self.bcy = bcy(1, self.j.start, self.j.stop, hlo)

    self.psi = (
      numpy.empty((
        ext(self.i, self.hlo, self.hlo).stop, 
        ext(self.j, self.hlo, self.hlo).stop
      ), real_t),
      numpy.empty((
        ext(self.i, self.hlo, self.hlo).stop, 
        ext(self.j, self.hlo, self.hlo).stop
      ), real_t)
    )

    self.C = (
      numpy.empty((
        ext(self.i, 1, 0).stop, 
        ext(self.j, self.hlo, self.hlo).stop
      ), real_t),
      numpy.empty((
        ext(self.i, self.hlo, self.hlo).stop, 
        ext(self.j, 1, 0).stop
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
    #pdb.set_trace()
    for t in range(nt):
      self.bcx.fill_halos(
        self.psi[self.n], ext(self.j, self.hlo, self.hlo).start, 
        ext(self.j, self.hlo, self.hlo).stop 
      )
      self.bcy.fill_halos(
        self.psi[self.n], ext(self.i, self.hlo, self.hlo).start, 
        ext(self.i, self.hlo, self.hlo).stop 
      )
      self.advop() 
      self.cycle()
  
#listing08
class Cyclic(object):
  # ctor
  def __init__(self, d, i_start, i_stop, hlo):
    i  = slice(i_start, i_stop)
    self.d = d
    self.left_halo = slice(i.start-hlo, i.start    )
    self.rght_edge = slice(i.stop -hlo, i.stop     )
    self.rght_halo = slice(i.stop,      i.stop +hlo)
    self.left_edge = slice(i.start,     i.start+hlo)

  # method invoked by the solver
  def fill_halos(self, psi, j_start, j_stop):
    j = slice(j_start, j_stop) 
    psi[pi(self.d, self.left_halo, j)] = (
      psi[pi(self.d, self.rght_edge, j)]
    )
    psi[pi(self.d, self.rght_halo, j)] = (
      psi[pi(self.d, self.left_edge, j)]
    )

#listing09
@autojit
def f(psi_l, psi_r, C):
  #pdb.set_trace()
  return (
    (C + abs(C)) * psi_l + 
    (C - abs(C)) * psi_r
  ) / 2

#listing10
@autojit
def donorcell(d, psi, C, i, j):
  #pdb.set_trace()
  return (
    f(
      psi[pi(d, i,     j)], 
      psi[pi(d, sl_add(i,1), j)], 
        C[pi(d, sl_add(i,0), j)]
    ) - 
    f(
      psi[pi(d, sl_add(i,-1), j)], 
      psi[pi(d, i,     j)], 
        C[pi(d, sl_add(i,-1), j)]
    ) 
  )

#listing11
@autojit
def donorcell_op(psi, n, C, i, j):
  #pdb.set_trace()
  psi[n+1][i,j] = (psi[n][i,j] 
    - donorcell(0, psi[n], C[0], i, j)
    - donorcell(1, psi[n], C[1], j, i)
  )

#listing12
class Donorcell(Solver):
  def __init__(self, bcx, bcy, nx, ny):
    Solver.__init__(self, bcx, bcy, nx, ny, 1)

  def advop(self):
    donorcell_op(
      self.psi, self.n, 
      self.C, self.i.start, self.i.stop, self.j.start, self.j.stop
    )

#listing13
#@autojit
def frac(nom, den):
  return numpy.where(den > 0, nom/den, 0)

#listing14
@autojit
def a_op(d, psi, i, j):
  #pdb.set_trace()
  return frac(
    psi[pi(d, sl_add(i,1), j)] - psi[pi(d, i, j)],
    psi[pi(d, sl_add(i,1), j)] + psi[pi(d, i, j)]
  )

#listing15
@autojit
def b_op(d, psi, i, j):
  return frac( 
    psi[pi(d, sl_add(i,1), sl_add(j,1))] + psi[pi(d, i, sl_add(j,1))] -
    psi[pi(d, sl_add(i,1), sl_add(j,-1))] - psi[pi(d, i, sl_add(j,-1))],
    psi[pi(d, sl_add(i,1), sl_add(j,1))] + psi[pi(d, i, sl_add(j,1))] +
    psi[pi(d, sl_add(i,1), sl_add(j,-1))] + psi[pi(d, i, sl_add(j,-1))]
  ) / 2

#listing16
@autojit
def C_bar(d, C, i, j):
  return (
    C[pi(d, sl_add(i,1), sl_add(j,0))] + C[pi(d, i,  sl_add(j,0))] +
    C[pi(d, sl_add(i,1), sl_add(j,-1))] + C[pi(d, i,  sl_add(j,-1))] 
  ) / 4

#listing17
@autojit
def antidiff(d, psi, i, j, C):
  return (
    abs(C[d][pi(d, sl_add(i,0), j)]) 
    * (1 - abs(C[d][pi(d, sl_add(i,0), j)])) 
    * a_op(d, psi, i, j)
    - C[d][pi(d, sl_add(i,0), j)] 
    * C_bar(d, C[d-1], i, j)
    * b_op(d, psi, i, j)
  )
#listing18
class Mpdata(Solver):
  def __init__(self, n_iters, bcx, bcy, nx, ny):
    Solver.__init__(self, bcx, bcy, nx, ny, 1)
    self.im = slice(self.i.start-1, self.i.stop)
    self.jm = slice(self.j.start-1, self.j.stop)

    self.n_iters = n_iters
  
    self.tmp = [(
      numpy.empty(self.C[0].shape, real_t),
      numpy.empty(self.C[1].shape, real_t)
    )]
    if n_iters > 2:
      self.tmp.append((
        numpy.empty(self.C[0].shape, real_t),
        numpy.empty(self.C[1].shape, real_t)
      ))

  def advop(self):
    #pdb.set_trace()
    for step in range(self.n_iters):
      if step == 0:
        donorcell_op(
          self.psi, self.n, self.C, self.i, self.j
        )
      else:
        self.cycle()
        self.bcx.fill_halos(
          self.psi[self.n], ext(self.j, self.hlo, self.hlo).start, 
          ext(self.j, self.hlo, self.hlo).stop
        )
        self.bcy.fill_halos(
          self.psi[self.n], ext(self.i, self.hlo, self.hlo).start, 
          ext(self.i, self.hlo, self.hlo).stop 
        )
        if step == 1:
          C_unco, C_corr = self.C, self.tmp[0]
        elif step % 2:
          C_unco, C_corr = self.tmp[1], self.tmp[0]
        else:
          C_unco, C_corr = self.tmp[0], self.tmp[1]

        C_corr[0][sl_add(self.im,0), self.j] = (
          antidiff(0, self.psi[self.n], 
            self.im, self.j, C_unco)
          ) 
        self.bcy.fill_halos(C_corr[0], ext(self.i, 1, 0).start, ext(self.i, 1, 0).stop)
        
        C_corr[1][self.i, sl_add(self.jm,0)] = (
          antidiff(1, self.psi[self.n], 
            self.jm, self.i, C_unco)
          ) 
        self.bcx.fill_halos(C_corr[1], ext(self.j, 1, 0).start, ext(self.j, 1, 0).stop)

        donorcell_op(
          self.psi, self.n, C_corr, self.i, self.j
        )
#listing19
