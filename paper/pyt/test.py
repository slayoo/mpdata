import numpy as n

class h():
  pass
# TODO: operator overloading

def f(psi_l, psi_r, C):
  return n.where(C >= 0, psi_l * C, psi_r * C)

def donorcell_1D(psi, C, i, j):
  h = .5
  return (
    f(psi[i,   j], psi[i+1, j], C[i+h, j]) - 
    f(psi[i-1, j], psi[i,   j], C[i-h, j])
  )

def donornell_2D(psi, n, C, i, j):
  psi[n+1][i,j] = (psi[n][i,j]
    - donorcell_1D(psi[n], c[0], i, j)
    - donorcell_1D(psi[n].swap, c[1], j, i)
  )

class Mpdata(object):
  def __init__(self, n_iters):
    self.n_steps = n_iters
    self.n_halos = n_iters


class Solver_2D(object):
  def __init__(self, adv, nx, ny):
    self.n = 0
    self.psi = []
    for 
      self.psi.append(n.zeros((nx, ny))) #TODO: uninitialised

  def solve(self, nt):
    #print self.psi
    pass

  def state(self):
    return self.psi[self.n]




nx, ny = 32, 32

slv1 = Solver_2D(Mpdata(1), nx, ny)
slv1.solve(20)

print slv1.state()

slv2 = Solver_2D(Mpdata(2), nx, ny)
slv2.solve(20)
