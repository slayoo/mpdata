#listing00
# code licensed under the terms of GNU GPL v3
# copyright holder: University of Warsaw
#listing01
real_t = 'float64'
#listing02
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

    # u nas bcx jest zawsze cycling, uzupelnia halo w psi
    self.bcx = bcx(0, self.i, hlo)
    self.bcy = bcy(1, self.j, hlo)

    self.psi = (
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t),
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t) 
    )
    #chyba mogloby byc w pierwszym (nx+1,ny) i analogicznie, tylko wtedy w donorcel trzeba by inaczej j dla C zapisywac (musialoby byc j-1)
    #moglby byc (nx+1,..) jesli hlf byloby Shift(0,1)
    self.C = (
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t),
      numpy.empty((nx+2*hlo,   ny+2*hlo), real_t)
    )

  # accessor methods
  def state(self):
    return self.psi[self.n][self.i, self.j]
  def Cx(self):
    return self.C[0][self.i, self.j]
  def Cy(self):
    return self.C[1][self.i, self.j]

   # integration logic
  def solve(self, nt):
    self.bcx.fill_halos_vel(self.C)
    self.bcy.fill_halos_vel(self.C)
    for t in range(nt):
      for s in range(self.adv.n_steps):
        self.bcx.fill_halos(self.psi[self.n])
        self.bcy.fill_halos(self.psi[self.n])
        self.adv.op_2D(self.psi, self.n, 
          self.C, self.i, self.j, s, 
           self.bcx, self.bcy)
        #nie rozumiem!
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

  def fill_halos_vel(self, C):
    for i in range(len(C)):
      self.fill_halos(C[i])
    
    
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
    #mim zdaniem wystarczy zawsze jedno halo
    self.n_halos = 1
    hlo = self.n_halos
    #czy w konstruktorze powinno to byc?
    self.tmp0 = (
      numpy.empty((nx+2*hlo, ny+2*hlo), real_t),
      numpy.empty((nx+2*hlo,  ny+2*hlo), real_t)
    )
    if n_iters > 2:
      self.tmp1 = (
        numpy.empty(( nx+2*hlo, ny+2*hlo), real_t),
        numpy.empty(( nx+2*hlo,   ny+2*hlo), real_t)
        )

  def op_2D(self, psi, n, C, i, j, step, bcx, bcy):
    if step == 0:
      donorcell_2D(psi, n, C, i, j)
    else:
      if step == 1:
        C_unco, C_corr = C, self.tmp0
      elif step % 2:
        C_unco, C_corr = self.tmp1, self.tmp0
      else:
        C_unco, C_corr = self.tmp0, self.tmp1

#jesli hlh = Shift(0,1), to by to chyba nie bylo potrzebne
      im = i - one # przeniesc do konstruktora
      jm = j - one

      C_corr[0][im+hlf, j] = (
        antidiff_2D(0, psi[n], im, j, C_unco)) 
      C_corr[1][i, jm+hlf] = (
        antidiff_2D(1, psi[n], jm, i, C_unco)) 

      bcx.fill_halos_vel(C_corr)
      bcy.fill_halos_vel(C_corr)

      donorcell_2D(psi, n, C_corr, i, j)


      
