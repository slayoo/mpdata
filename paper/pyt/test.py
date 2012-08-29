# -- listing 2.0
import numpy #TODO: as np?

# -- listing 2.3

class half():
  # TODO? only if isinstance(slice)?
  def __radd__(self, a):
    return slice(a.start + 1, a.stop + 1)
  def __rsub__(self, a):
    return a

class one():
  # TODO? only if isinstance(slice)?
  def __radd__(self, a):
    return slice(a.start + 1, a.stop + 1)
  def __rsub__(self, a):
    return slice(a.start - 1, a.stop - 1)

h = half()
l = one() # TODO: is it possible not to use this hack???

# -- listing 2.6

def f(psi_l, psi_r, C):
  return numpy.where(C >= 0, psi_l * C, psi_r * C)

# -- listing 2.7

def donorcell_1D(psi, C, i, j):
  return (
    f(psi[i,   j], psi[i+l, j], C[i+h, j]) - 
    f(psi[i-l, j], psi[i,   j], C[i-h, j])
  )

# -- listing 2.8

def donorcell_2D(psi, n, C, i, j):
  psi[n+1][i,j] = (psi[n][i,j]
    - donorcell_1D(psi[n][:], C[0][:], i, j)
    - donorcell_1D(psi[n][:], C[1][:], i, j)
  )

# -- listing 2.13

class Mpdata(object):
  def __init__(self, n_iters):
    self.n_steps = n_iters
    self.n_halos = n_iters

  def op_2D(self, psi, n, C, i, j, step):
    donorcell_2D(psi, n, C, i, j)

# -- listing 2.15

class Solver_2D(object):
  def __init__(self, adv, nx, ny):
    self.adv = adv
    self.n = 0
    self.psi = []
    self.C = []
    self.i = slice(adv.n_halos, nx + adv.n_halos)
    self.j = slice(adv.n_halos, ny + adv.n_halos)

    hlo = adv.n_halos
    for i in range(2):
      self.psi.append(numpy.empty((nx + 2 * hlo, ny + 2 * hlo))) 
    self.C.append(numpy.empty((nx + 1 + 2 * hlo, ny + 2 * hlo)))
    self.C.append(numpy.empty((nx + 2 * hlo, ny + 1 + 2 * hlo)))

  def solve(self, nt):
    for t in range(nt):
      for s in range(self.adv.n_steps):
        # TODO: fill halos
        print "t=%d s=%d" % (t, s)
        self.adv.op_2D(self.psi, self.n, self.C, self.i, self.j, s)
#        self.n  = (self.n + 1) % 2

  def state(self):
    return self.psi[self.n]
  def Cx(self):
    return self.C[0]
  def Cy(self):
    return self.C[1]

# -- listing 2.16

nx, ny = 10, 5

slv1 = Solver_2D(Mpdata(1), nx, ny)

slv1.state()[:] = 0
slv1.state()[0:2,0:3] = 1
slv1.Cx()[:] = .5
slv1.Cy()[:] = .5

print slv1.state()

nt = 3

slv1.solve(nt)
print slv1.state()

slv2 = Solver_2D(Mpdata(2), nx, ny)
slv2.solve(nt)
