import numpy
import time

Psi_l = numpy.array([[ 0.,  0.,  0.],
                     [ 0.,  1.,  0.],
                     [ 0.,  0.,  0.]])

Psi_r = numpy.array([[ 0.,  0.,  0.],
                     [ 1.,  0.,  0.],
                     [ 0.,  0.,  0.]])

CC = numpy.array([[ 1.,  1.,  1.],
                 [ 1.,  1.,  1.],
                 [ 1.,  1.,  1.]])

#Psi_l = numpy.array([0., 1., 0.])

#Psi_r = numpy.array([1., 0., 0.])

#CC = numpy.array([1., 1., 1.])

def f(psi_l, psi_r, C):
    return (
        (C + abs(C)) * psi_l +
        (C - abs(C)) * psi_r
        ) / 2

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

one = Shift(1,1)
sl = slice(2,5)
print sl+one

try:
    from numba import autojit, jit, double
    fast_f = autojit(f)
    #fast_f = jit(restype=double[:],argtypes=[double[:], double[:], double[:]])(f)
    t0 = time.time()
    output = fast_f(Psi_l, Psi_r, CC)
    print "time: numba", time.time() - t0
    print "output numba", output
except ImportError:
    print "Not using Numba nor Numbapro"

t1 = time.time()
output = f(Psi_l, Psi_r, CC)
print "time: numpy", time.time() - t1
print "output numpy", output


