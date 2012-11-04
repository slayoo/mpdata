#listing00
# code licensed under the terms of GNU GPL v3
# copyright holder: University of Warsaw
#listing01
real_t = 'float64'
#listing02
import numpy
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
hlf = Shift(1,0)
#listing05
class Solver_2D(object):
  # ctor
  def __init__(self, adv, n_iters, bcx, bcy, nx, ny):
    self.adv = adv(n_iters, nx, ny)
    self.n = 0

    hlo = self.adv.n_halos

    self.i = slice(hlo, nx + hlo)
    self.j = slice(hlo, ny + hlo)

    self.bcx = bcx(0, self.i, hlo)
    self.bcy = bcy(1, self.j, hlo)

    self.psi = (
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t),
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t) 
    )
    self.C = (
      numpy.empty((nx+1+2*hlo, ny+2*hlo), real_t),
      numpy.empty((nx+2*hlo,   ny+1+2*hlo), real_t)
    )

  # accessor methods
  def state(self):
    return self.psi[self.n][self.i, self.j]
  def Cx(self):
    return self.C[0]
  def Cy(self):
    return self.C[1]

  # integration logic
  def solve(self, nt):
    for t in range(nt):
      for s in range(self.adv.n_steps):
        self.bcx.fill_halos(self.psi[self.n])
        self.bcy.fill_halos(self.psi[self.n])
        self.adv.op_2D(self.psi, self.n, 
          self.C, self.i, self.j, s
        )
        self.n  = (self.n + 1) % 2 - 2
        
#listing06
def pi(d, *idx): 
  return (idx[d], idx[d-1])
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
def donorcell_2D(psi, n, C, i, j):
  psi[n+1][i,j] = (psi[n][i,j] 
    - donorcell(0, psi[n], C[0], i, j)
    - donorcell(1, psi[n], C[1], j, i)
  )
#listing11
def frac(nom, den):
  return numpy.where(den>0, nom/den, 0)
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
def antidiff_2D(d, psi, i, j, C):
  return (
    abs(C[d][pi(d, i+hlf, j)]) 
    * (1 - abs(C[d][pi(d, i+hlf, j)])) 
    * a_op(d, psi, i, j)
    - C[d][pi(d, i+hlf, j)] 
    * ( # zrobic f-cje vmean
      C[d-1][pi(d, i+one, j+hlf)] + 
      C[d-1][pi(d, i,     j+hlf)] +
      C[d-1][pi(d, i+one, j-hlf)] + 
      C[d-1][pi(d, i,     j-hlf)] 
    ) / 4
    * b_op(d, psi, i, j)
  )
#listing15
class Mpdata(object):
  def __init__(self, n_iters, nx, ny):
    self.n_steps = n_iters
    self.n_halos = n_iters
    hlo = self.n_halos
    self.tmp0 = (
      numpy.empty(( nx+1+2*hlo, ny+2*hlo), real_t),
      numpy.empty(( nx+2*hlo,   ny+1+2*hlo), real_t)
    )
    if n_iters > 2:
      self.tmp1 = (
        numpy.empty(( nx+1+2*hlo, ny+2*hlo), real_t),
        numpy.empty(( nx+2*hlo,   ny+1+2*hlo), real_t)
      )

  def op_2D(self, psi, n, C, i, j, step):
    if step == 0:
      donorcell_2D(psi, n, C, i, j)
    else:
      if step == 1:
        C_unco, C_corr = C, self.tmp0 #zmienic tmp0,1 na tuple tmp
      elif step % 2:
        C_unco, C_corr = self.tmp1, self.tmp0
      else:
        C_unco, C_corr = self.tmp0, self.tmp1
        
        
      im = slice(i.start - 1, i.stop) # przeniesc do konstruktora
      jm = slice(j.start - 1, j.stop)

      # not sure if this is ok, ask Sylwester
      # TO DO: should we write similar class to Shift??
      ih = slice(i.start - (self.n_steps - 1 - step), i.stop +  (self.n_steps - 1 - step)) 
      jh = slice(j.start -  (self.n_steps - 1 - step), j.stop +  (self.n_steps - 1 - step))

      print "rozmiar psi", psi[n].shape, im, jm, ih, jh, self.n_halos
      
      C_corr[0][im+hlf, jh] = (
        antidiff_2D(0, psi[n], im, jh, C_unco)) 
      C_corr[1][ih, jm+hlf] = (
        antidiff_2D(1, psi[n], jm, ih, C_unco)) 
      donorcell_2D(psi, n, C_corr, i, j)
#listing16
