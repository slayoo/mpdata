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

#    def left_halo TODO
#    def rght_halo TODO

    def rght_edge(self,n):
        return self.psi[n][-2*self.adv.halo : -self.adv.halo]
    
    def left_edge(self,n):
        return self.psi[n][self.adv.halo : 2*self.adv.halo]
             
    def fill_halos(self,n):
        self.psi[n][0 : self.adv.halo] = self.rght_edge(n)
        self.psi[n][-self.adv.halo :] = self.left_edge(n)

    def advect(self,adv,n,ns=0):
        adv.op(self.psi,self.v,n,(self.adv.halo, self.nxy + self.adv.halo),ns)

# derived class implementing the upstream algorithm
class AdvUpstream(object):

    num_steps = 1
    stencil_extent = 3
    time_levels = 2
    halo = (stencil_extent-1)/2

    def f(self,psi_l,psi_r,U):
        return N.where(U >=0, psi_l * U, psi_r * U)

    def op(self,psi_tab,v,n,i_range,ns):
        i_p, i_k = i_range
        psi_tab[n+1][i_p:i_k] -= (self.f(psi_tab[n][i_p:i_k],psi_tab[n][i_p+1:i_k+1],v[i_p+1:i_k+1])
                               - self.f(psi_tab[n][i_p-1:i_k-1],psi_tab[n][i_p:i_k],v[i_p:i_k]))

nf = NetCDFFile(sys.argv[1], 'r+')
#wczytywanie atrybutow z pliku netcdf
nx = int(nf.nx)
ny = int(nf.ny)
no = int(nf.no)
np = nf.np
nt = nf.nt
cx = nf.Cx
cy = nf.Cy
#wczytywania tablicy skalara psi
nv = nf.variables['psi'][:]
  
#zakladam chwilowo, ze adwekcja jest tylko upsream
advop = AdvUpstream()
    
psi = []
# jesli mialoby to byc rownolegle, to nx -> ix_max-ix_min (chyba)
for n in range(advop.time_levels):
    psi.append(N.zeros((nx+2*advop.halo,ny+2*advop.halo)))
#tworzenie kopii odwroconej
psi_sw = []
for n in range(advop.time_levels):
    psi_sw.append(psi[n].swapaxes(0,1))

#pole preskosci
v = []
#w wyzszych iteracjach potrzbne halo
for n in range(2): #potrzebuje predkosc w x i y
    v.append(N.zeros((nx+2*advop.halo,ny+2*advop.halo))) # nie potrzeba tak naprawde dodawac halo. rozmiar wewnetrznej tablicy musi by tylko o 1 wiekszy niz w psi, ale tak mi wygodniej chwilowo
        
#tworzenie kopii odwroconej
v_sw = []
for n in range(2):
    v_sw.append(v[n].swapaxes(0,1))

#war pocztkowy przeniesc?
psi[0][advop.halo : nx+advop.halo, advop.halo : ny+advop.halo] = nv
v[0][1:nx+2,1:ny+2] = cx
v[1][1:nx+2,1:ny+2] = cy

#tworze obiekt osobny dla kierunku X i Y? Nie wiem, czy tak najlepiej.??
advop_X = Adv(advop,psi,v[0])
advop_Y = Adv(advop,psi_sw,v_sw[1])

for t in range(1,nt+1):
    n = 0 #aby pisac psi[n+1] = psi[n]+...
    #wypelnianie hala w dwoch wymiarach
    advop_X.fill_halos(n)
    advop_Y.fill_halos(n)
    psi[n+1][:] = psi[n]
    advop_X.advect(advop,n)
    advop_Y.advect(advop,n)
    #print "t, psi po adwekcji", t, psi
    psi[n][:] = psi[n+1]
        
    if t%no == 0:
         nf.variables['psi'][t/no,:] = psi[0][advop.halo : nx+advop.halo, advop.halo : ny+advop.halo].astype('f')
nf.close()
