# array handling: Numpy
import numpy as N

# I/O: netCDF-4 (via Scientific Python)
from Scientific.IO.NetCDF import NetCDFFile

# for handling command line arguments
import sys

# abstract class defining the advection operator 
class Adv(object):
    def __init__(self,adv,psi_tab,v_tab):
        self.adv = adv
        self.psi = psi_tab
        self.v = v_tab
        self.nxy = self.psi[0].shape[0] - 2*self.adv.halo

    def rght_edge(self,n):
        return self.psi[n][-2*self.adv.halo:-self.adv.halo]
    
    def left_edge(self,n):
        return self.psi[n][self.adv.halo:2*self.adv.halo]
             
    def fill_halos(self,n):
        self.psi[n][0:self.adv.halo] = self.rght_edge(n)
        self.psi[n][-self.adv.halo:] = self.left_edge(n)

    def advect(self,adv,n,ns=0):
        adv.op(self.psi,self.v,n,(self.adv.halo, self.nxy + self.adv.halo),ns)

# derived class implementing the upstream algorithm
class AdvUpstream(object):
    num_steps = 1
    stencil_extent = 3
    time_levels = 2
    halo = (stencil_extent-1)/2

    # e.g. eq. (3 a-d) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480)
    def f(self,psi_l,psi_r,U):
        return N.where(U >=0, psi_l * U, psi_r * U)

    # e.g. eq. (2) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480)
    def op(self,psi_tab,v,n,i_range,ns):
        i_p, i_k = i_range
        psi_tab[n+1][i_p:i_k] -= (
          self.f(psi_tab[n][i_p  : i_k  ],psi_tab[n][i_p+1 : i_k+1], v[i_p+1: i_k+1]) -
          self.f(psi_tab[n][i_p-1: i_k-1],psi_tab[n][i_p   : i_k  ], v[i_p  : i_k  ])
        )

# handling command line argument
if (len(sys.argv) != 2) : 
  raise Exception('expecting 1 argument - a netCDF file name')

# opening and reading-in data from the netCDF file
nf = NetCDFFile(sys.argv[1], 'r+')
nx = int(nf.nx)
ny = int(nf.ny)
no = int(nf.no)
np = nf.np
nt = nf.nt
  
# instantiationg the advection operator
advop = AdvUpstream()
    
# allocating memory for psi (at two time levels)
psi = []
for n in range(advop.time_levels):
    psi.append(N.zeros((nx+2*advop.halo, ny+2*advop.halo)))

# allocating memory for the velocity field (x & y)
vel = []
for n in range(2): 
    vel.append(N.zeros((nx+2*advop.halo, ny+2*advop.halo))) 
        
# helper views for dimension-independent logic
psi_sw = []
for n in range(advop.time_levels):
    psi_sw.append(psi[n].swapaxes(0,1))
vel_sw = []
for n in range(2):
    vel_sw.append(vel[n].swapaxes(0,1))

# filling psi and vel with data from netCDF
psi[0][advop.halo : nx+advop.halo, advop.halo : ny+advop.halo] = (
  nf.variables['psi'][:]
)
vel[0][1:nx+2,1:ny+2] = nf.Cx
vel[1][1:nx+2,1:ny+2] = nf.Cy

# dimension-independency logic contd.
advop_X = Adv(advop,psi,vel[0])
advop_Y = Adv(advop,psi_sw,vel_sw[1])

# integration loop
for t in range(1,nt+1):
    n = 0 # for human readibility :)

    # filling the halos
    advop_X.fill_halos(n)
    advop_Y.fill_halos(n)

    # advecting in each dimension
    psi[n+1][:] = psi[n]
    advop_X.advect(advop,n)
    advop_Y.advect(advop,n)

    # outputtting to the netCDF
    if t%no == 0:
        nf.variables['psi'][t/no,:] = (
          psi[n+1][advop.halo:nx+advop.halo, advop.halo:ny+advop.halo].astype('f')
        )

    # cycling the pointers
    psi[n][:] = psi[n+1]

# closing the netCDF file
nf.close()
