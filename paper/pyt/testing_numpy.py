"""Unit test for listing.py using numpy function for testing"""
from listings import *
import numpy

knownValue = [[3,3,1,0.5,0.5,
               numpy.array([[0., 0, 0],
                            [0,  1, 0],
                            [0,  0, 0]]),
               numpy.array([[0.,  0,   0],
                            [0,   0.,  0.5],
                            [0,   0.5, 0]]),
               numpy.array([[0, 0,   0],
                            [0, 0.,  0.5],
                            [0, 0.5, 0]])],
              [3,3,1,0.2,0.2,
               numpy.array([[0., 0, 0],
                            [0,  1, 0],
                            [0,  0, 0]]),
               numpy.array([[0, 0,   0],
                            [0, 0.6, 0.2],
                            [0, 0.2, 0.]]),
               numpy.array([[0, 0,     0],
                            [0, 0.64,  0.18],
                            [0, 0.18,  0]])],
              [3,3,1,0.1,0.5,
               numpy.array([[0., 0, 0],
                            [0,  1, 0],
                            [0,  0, 0]]),
               numpy.array([[0, 0,   0],
                            [0, 0.4, 0.5],
                            [0, 0.1, 0]]),
               numpy.array([[0, 0,       0],
                            [0, 0.4068,  0.5011],
                            [0, 0.0921,  0]])]
              ]

def testMpdata(nx, ny, nt, cx, cy, n_iters, psi_in, psi_out):
    slv = Solver_2D(Mpdata, n_iters, Cyclic, Cyclic, nx, ny)
    slv.state()[:] = psi_in
    slv.Cx()[:] = cx
    slv.Cy()[:] = cy
    slv.solve(nt)
    print "test dla parametrow: cx, cy, nt, n_iters", cx, cy, nt, n_iters
    print "psi rozw", slv.state(), slv.state() - psi_out
    print "test", numpy.allclose(slv.state(), psi_out, atol=1.e-5), '\n'

def main():
    for nx, ny, nt, cx, cy, psi_in, psi_out1, psi_out2 in knownValue:
        testMpdata(nx, ny, nt, cx, cy, 1, psi_in, psi_out1)
        testMpdata(nx, ny, nt, cx, cy, 2, psi_in, psi_out2)

main()
