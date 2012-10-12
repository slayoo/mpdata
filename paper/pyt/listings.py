#listing00
import numpy 
#listing01
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
#listing02
one = Shift(1,1) 
hlf = Shift(1,0)
#listing02
def pi(d, idx): 
  return (idx[d], idx[d-1])
#listing03
class Solver_2D(object):
  def __init__(self, adv, bcx, bcy, nx, ny):
    self.adv = adv
    self.n = 0

    hlo = adv.n_halos

    self.i = slice(hlo, nx + hlo)
    self.j = slice(hlo, ny + hlo)

    self.bcx = bcx(0, self.i, hlo)
    self.bcy = bcy(1, self.j, hlo)

    self.psi = (
      numpy.empty((nx + 2 * hlo, ny + 2 * hlo)),
      numpy.empty((nx + 2 * hlo, ny + 2 * hlo)) 
    )
    self.C = (
      numpy.empty((nx + 1 + 2 * hlo, ny + 2 * hlo)),
      numpy.empty((nx + 2 * hlo, ny + 1 + 2 * hlo))
    )

  def solve(self, nt):
    for t in range(nt):
      for s in range(self.adv.n_steps):
        print "t=%d s=%d" % (t, s)
        self.bcx.fill_halos(self.psi[self.n])
        self.bcy.fill_halos(self.psi[self.n])
        self.adv.op_2D(self.psi, self.n, self.C, self.i, self.j, s)
        self.n  = (self.n + 1) % 2 - 2

  def state(self):
    return self.psi[self.n][self.i,self.j]
  def Cx(self):
    return self.C[0]
  def Cy(self):
    return self.C[1]
#listing04
class Cyclic(object):
  def __init__(self, d, i, hlo): #TODO: i & j not needed in principle...
    self.d = d
    self.left_halo = slice(i.start - hlo, i.start)
    self.rght_edge = slice(i.stop - hlo, i.stop)
    self.rght_halo = slice(i.stop, i.stop + hlo)
    self.left_edge = slice(i.start, i.start + hlo)

  def fill_halos(self, psi):
    psi[pi(self.d, [self.left_halo, slice(None)])] = psi[pi(self.d, [self.rght_edge, slice(None)])]
    psi[pi(self.d, [self.rght_halo, slice(None)])] = psi[pi(self.d, [self.left_edge, slice(None)])]
#listing05
def f(psi_l, psi_r, C): 
  return numpy.where(C >= 0, psi_l * C, psi_r * C)
#listing06
def donorcell_1D(d, psi, C, i, j): 
  return (
    f(psi[pi(d, [i,     j])], psi[pi(d, [i+one, j])], C[pi(d, [i+hlf, j])]) - 
    f(psi[pi(d, [i-one, j])], psi[pi(d, [i,     j])], C[pi(d, [i-hlf, j])]) 
  )
#listing07
def donorcell_2D(psi, n, C, i, j):
  psi[n+1][i,j] = (psi[n][i,j] 
    - donorcell_1D(0, psi[n], C[0], i, j)
    - donorcell_1D(1, psi[n], C[1], j, i)
  )
#listing08
class Mpdata(object):
  def __init__(self, n_iters):
    self.n_steps = n_iters
    self.n_halos = n_iters
  	

  def op_2D(self, psi, n, C, i, j, step):
    donorcell_2D(psi, n, C, i, j)
#listing09
