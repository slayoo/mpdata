// garbage collection
#include <boost/ptr_container/ptr_vector.hpp>
using boost::ptr_vector;

// working in 2 dimensions, with single precision
const int ndim = 2;
typedef float real_t;

// OpenMP for shared-memory parallelisation
#if defined(_OPENMP)
#  include <omp.h>
#  define BZ_THREADSAFE
#endif

// Blitz++ for arrays and array expressions
#include <blitz/array.h>
typedef blitz::Array<real_t, ndim> arr_t;
typedef ptr_vector<arr_t> vec_arr_t;
typedef blitz::Range rng_t;
typedef blitz::RectDomain<ndim> idx_t;
typedef blitz::TinyVector<rng_t, ndim> tv_t;

// netCDF-4 for input and output
#include <netcdf>
typedef long unsigned int siz_t;
typedef std::vector<siz_t> pos_t;

// some helper constructs
struct idx_ij : idx_t { idx_ij(const rng_t &i, const rng_t &j) : idx_t(tv_t(i,j)) {} }; 
struct idx_ji : idx_t { idx_ji(const rng_t &j, const rng_t &i) : idx_t(tv_t(i,j)) {} }; 
#define decltype_return(expr) -> decltype(expr) { return expr; } // requires C++11
#define error(msg) { std::cerr << msg << std::endl; throw std::exception(); }

// abstract class defining the advection operator functor
class adv {
  public: static const int ntlev = 2; // two-time leve schemes
  public: virtual void operator()(const vec_arr_t &psi, const arr_t &C, int n, int s) = 0;
  public: virtual rng_t rng_psi(int n) = 0;
  public: virtual rng_t rng_vel(int n) = 0;
  public: virtual idx_t left_halo() = 0;
  public: virtual idx_t rght_halo() = 0;
  public: virtual idx_t left_edge() = 0;
  public: virtual idx_t rght_edge() = 0;
  public: virtual int n_steps() { return 1; }
  protected: static const int phalf = 1, mhalf = 0; // for Arakawa-C staggered grid
};

// derived abstract class with the 1D->2D logic
template <class idx> class adv_idx : public adv {
  private: rng_t rng(int n, int p, int np) { return rng_t(p*n/np, (p+1)*n/np-1); }
  public: adv_idx(int nx, int ny, int px, int npx, int py, int npy, int halo) :
    iall(0,nx-1), jall(0,ny-1), i(rng(nx, px, npx)), j(rng(ny, py, npy)), halo(halo) {}
  private: const int halo;
  public: rng_t rng_psi(int n) { return rng_t(0 - halo, n - 1 + halo); }
  public: rng_t rng_vel(int n) { return rng_t(0 - halo - mhalf, n - 1 + halo + phalf); }
  public: idx_t left_halo() { return idx(rng_t(iall.first() - halo, iall.first() - 1), jall); }
  public: idx_t rght_halo() { return idx(rng_t(iall.last() + 1, iall.last() + halo), jall); }
  public: idx_t left_edge() { return idx(rng_t(iall.first(), iall.first() + halo - 1), jall); }
  public: idx_t rght_edge() { return idx(rng_t(iall.last() - halo + 1, iall.last()), jall); }
  protected: const rng_t i, j; // TODO: czy i,j jako argumenty do op() nie upro¶ci³yby?...
  private: const rng_t iall, jall;
};

// derived class implementing the upstream aka upwind aka donor-cell algorithm 
template <class idx> class adv_upstream : public adv_idx<idx> {
  public: adv_upstream(int nx, int ny, int px, int npx, int py, int npy) : 
    adv_idx<idx>(nx, ny, px, npx, py, npy, 1) {}

  /// eq. (3 a-d) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480)
  protected: template <class arg_t> static auto F(arg_t psi_l, arg_t psi_r, arg_t U)
    decltype_return(max(0,U) * psi_l + min(0,U) * psi_r)

  /// eq. (2) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) 
  public: void operator()(const vec_arr_t &psi, const arr_t &vel, int n, int s) {
    const rng_t &i = this->i, &j = this->j;
    psi[n+1](idx(i,j)) -= (
      F(psi[n](idx(i,  j)), psi[n](idx(i+1,j)), vel(idx(i + adv::phalf,j))) -
      F(psi[n](idx(i-1,j)), psi[n](idx(i,  j)), vel(idx(i - adv::mhalf,j)))
    );
  }
};

// derived class implementing the MPDATA algorithm 
/*
template <class idx> class adv_mpdata : public adv_idx<idx> {
  private: const adv_upstream<idx> upstream;
  public: adv_mpdata(int nx, int ny, int px, int npx, int py, int npy) :
    adv_idx<idx>(nx, ny, px, npx, py, npy, TODO halo(iter)) {}
    upstream(nx, ny, px, npx, py, npy)
  {}

// TODO: vel pointer
// TODO: halo corners!
// TODO: reference
  public: int n_steps() { return 2; }
  public: void operator()(const vec_arr_t &psi, const arr_t &vel, int n, int s) {
    switch (s) {
      default:
        // TODO
      case 1: 
        // TODO
      case 0: 
        return upstream(psi, vel, n, s);
    }
  }
};
*/

// execution starts here
int main(int ac, char* av[]) {

  // reading in simulation parameters from the netCDF file
  if (ac != 2) error(av[0] << " expects one argument - a netCDF file name")
  siz_t nt, nx, ny, np, no;
  real_t Cx, Cy;
  netCDF::NcFile nf(std::string(av[1]), netCDF::NcFile::write);
  netCDF::NcVar nv = nf.getVar("psi");
  nf.getAtt("nx").getValues(&nx);
  nf.getAtt("ny").getValues(&ny);
  nf.getAtt("np").getValues(&np);
  nf.getAtt("nt").getValues(&nt);
  nf.getAtt("no").getValues(&no);
  nf.getAtt("Cx").getValues(&Cx);
  nf.getAtt("Cy").getValues(&Cy);

  // chosing the number of threads to use
#if defined(_OPENMP)
  omp_set_num_threads(std::min(np, nx)); // np > nx does not make sense
#endif

  // memory allocation: advection functors
  std::vector<ptr_vector<adv>> advop(np); 
  for (int p = 0; p < np; ++p) { 
    advop[p].push_back(new adv_upstream<idx_ij>(nx, ny, p, np, 0, 1)); // x
    advop[p].push_back(new adv_upstream<idx_ji>(ny, nx, 0, 1, p, np)); // y
  }

  // memory allocation: scalar fields at two time levels (n, n+1)
  vec_arr_t psi(ndim); 
  for (int i=0; i<adv::ntlev; ++i) 
    psi.push_back(new arr_t(advop[0][0].rng_psi(nx), advop[0][0].rng_psi(ny)));

  // memory allocation: Courant number fields for two dimensions (x,y)
  vec_arr_t vel(ndim); 
  for (int i=0; i<ndim; ++i) 
    vel.push_back(new arr_t(advop[0][0].rng_vel(nx), advop[0][0].rng_vel(ny)));
  vel[0] = Cx;
  vel[1] = Cy;

  // reading the initial condition (due to presence of haloes the data to be stored 
  // is not contiguous, hence looping over one dimension) dims = {t,x,y}
  for (
    pos_t ncdf_off({0,0,0}), ncdf_cnt({1,1,size_t(ny)}); 
    ncdf_off[1] < nx; 
    ++ncdf_off[1]
  ) nv.getVar(ncdf_off, ncdf_cnt, psi[0](ncdf_off[1], rng_t(0,ny)).dataFirst());

  // integration
  for (siz_t t = 1; t <= nt; ++t) {
    const int n = 0; // that's just to resemble the equations

    for (int s = 0; s < advop[0][0].n_steps(); ++s) {
      // filling halos 
      for (int d = 0; d < ndim; ++d) {
        psi[n](advop[0][d].left_halo()) = psi[n](advop[0][d].rght_edge());
        psi[n](advop[0][d].rght_halo()) = psi[n](advop[0][d].left_edge());
      }

      // advecting in each dimension sharing the workload among threads
      psi[n+1] = psi[n];
#pragma omp parallel for
      for (int p = 0; p < np; ++p) // threads
        for (int d = 0; d < ndim; ++d) // dimensions
          advop[p][d](psi, vel[d], n, s);
    }

    // outputting
    if ((t % no) == 0) for (
      pos_t ncdf_off({t/no,0,0}), ncdf_cnt({1,1,ny}); 
      ncdf_off[1] < nx; 
      ++ncdf_off[1]
    ) nv.putVar(ncdf_off, ncdf_cnt, psi[n+1](ncdf_off[1], rng_t(0,ny)).dataFirst());

    // cycling the pointers to the arrays
    blitz::cycleArrays(psi[n], psi[n+1]);
  }  
}
