// garbage collection
#include <boost/ptr_container/ptr_vector.hpp>
using boost::ptr_vector;

// working in 2 dimensions, with double precision
const int ndim = 2;
typedef double real_t;

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
typedef size_t siz_t;
typedef std::vector<siz_t> pos_t;

// some helper constructs
struct idx_ij : idx_t { idx_ij(const rng_t &i, const rng_t &j) : idx_t(tv_t(i,j)) {} }; 
struct idx_ji : idx_t { idx_ji(const rng_t &j, const rng_t &i) : idx_t(tv_t(i,j)) {} }; 
#define decltype_return(expr) -> decltype(expr) { return expr; } // requires C++11
#define error(msg) { std::cerr << msg << std::endl; throw std::exception(); }

// abstract class defining the advection operator functor
class adv {
  public: static const int ntlev = 2; // two-time leve schemes
  public: virtual void operator()(
    const vec_arr_t &psi, const arr_t &C, 
    const rng_t &i, const rng_t &j, const int n, const int s
  ) = 0;
  public: virtual rng_t rng_psi() = 0;
  public: virtual rng_t rng_vel() = 0;
  public: virtual idx_t left_halo() = 0;
  public: virtual idx_t rght_halo() = 0;
  public: virtual idx_t left_edge() = 0;
  public: virtual idx_t rght_edge() = 0;
  public: virtual int n_steps() { return 1; }
  protected: static const int ph = 1, mh = 0; // for Arakawa-C staggered grid
};

// derived abstract class with the 1D->2D logic
template <class idx> class adv_idx : public adv 
{
  private: const int nx, ny, halo;

  public: adv_idx(const int nx, const int ny, const int halo) : nx(nx), ny(ny), halo(halo) {}

  public: rng_t rng_psi() 
  { return rng_t(0 - halo, nx - 1 + halo); }

  public: rng_t rng_vel() 
  { return rng_t(0 - halo - mh, nx - 1 + halo + ph); }

  public: idx_t left_halo() 
  { return idx(rng_t(- halo, -1), rng_t(-halo, ny - 1 + halo)); }

  public: idx_t rght_halo() 
  { return idx(rng_t(nx, nx - 1 + halo), rng_t(-halo, ny - 1 + halo)); }

  public: idx_t left_edge() 
  { return idx(rng_t(0, halo - 1), rng_t(-halo, ny - 1 + halo)); }

  public: idx_t rght_edge() 
  { return idx(rng_t(nx - halo, nx - 1), rng_t(-halo, ny - 1 + halo)); }
};

// derived class implementing the upstream aka upwind aka donor-cell algorithm 
template <class idx> class adv_upstream : public adv_idx<idx> 
{
  public: adv_upstream(const int nx, const int ny) : adv_idx<idx>(nx, ny, 1) {}

  // eq. (3 a-d) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480)
  protected: template <class a1_t, class a2_t, class a3_t> static auto F(a1_t psi_l, a2_t psi_r, a3_t U)
    decltype_return(max(0,U) * psi_l + min(0,U) * psi_r)

  // eq. (2) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) 
  public: void operator()(
    const vec_arr_t &psi, const arr_t &vel, 
    const rng_t &i, const rng_t &j, const int n, const int
  ) {
    psi[n+1](idx(i,j)) -= (
      F(psi[n](idx(i,  j)), psi[n](idx(i+1,j)), vel(idx(i + adv::ph,j))) -
      F(psi[n](idx(i-1,j)), psi[n](idx(i,  j)), vel(idx(i - adv::mh,j)))
    );
  }
};

template <class idx> class adv_mpdata : public adv_upstream<idx> 
{
  public: adv_mpdata(int nx, int ny) : adv_upstream<idx>(nx, ny) {}
  public: int n_steps() { return 2; }

  // preventing zeros from entering denominators
  protected: template <class a1_t, class a2_t> static auto frac(a1_t num, a2_t den)
    decltype_return(where(den > real_t(0), num / den, real_t(0)))

  // eq. (17a) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) 
  protected: static auto A(const arr_t &psi, const rng_t &i, const rng_t &j)
    decltype_return(frac(
      psi(i+1,j) - psi(i,j),
      psi(i+1,j) + psi(i,j)
    ))

  // eq. (17b) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) 
  protected: static auto B(const arr_t &psi, const rng_t &i, const rng_t &j)
    decltype_return(real_t(.5) * frac(
      psi(i+1,j+1) + psi(i,j+1) - psi(i+1,j-1) - psi(i,j-1),
      psi(i+1,j+1) + psi(i,j+1) + psi(i+1,j-1) + psi(i,j-1)
    ))

  // eq. (29b) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480)
  protected: static auto vmean(const arr_t &V, const rng_t &i, const rng_t &j)
    decltype_return(real_t(.25) * (
      V(i + 1,j + adv::ph) +
      V(i + 0,j + adv::ph) +
      V(i + 1,j - adv::ph) +
      V(i + 0,j - adv::ph)
    ))

  // eq. (29a) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) (for 2D case)
  protected: static auto U(const arr_t &psi, const arr_t &V, const rng_t &i, const rng_t &j)
    decltype_return(
      (abs(V(i + adv::ph, j)) - .5 * pow(V(i + adv::ph, j),2)) * A(psi, i, j)
      - vmean(V, i, j) * V(i + adv::ph, j) * B(psi, i, j)
    )

  public: void operator()(
    const vec_arr_t &psi, const arr_t &vel, 
    const rng_t &i, const rng_t &j, const int n, const int step
  ) {
    switch (step) {
      case 0: return adv_upstream<idx>::operator()(psi, vel, i, j, n, step);
      case 1: psi[n+1](idx(i,j)) -= (
        this->F(psi[n](idx(i,  j)), psi[n](idx(i+1,j)), U(psi[n], vel, i, j)) -
        this->F(psi[n](idx(i-1,j)), psi[n](idx(i,  j)), U(psi[n], vel, i-1, j))
      );
      return;
      default: assert(false);
    }
  }
};

// execution starts here
int main(int ac, char* av[]) 
{
  const int n = 0, x = 0, y = 1; // for human readibility :)

  // reading in simulation parameters from the netCDF file
  if (ac != 2) error(av[0] << " expects one argument - a netCDF file name")
  int nt, nx, ny, np, no;
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
  ptr_vector<adv> advop(ndim); 
  bool mpdata = true; // TODO: opcja linii komend?
  if (mpdata) {
    advop.push_back(new adv_mpdata<idx_ij>(nx, ny)); // x 
    advop.push_back(new adv_mpdata<idx_ji>(ny, nx)); // y 
  } else {
    advop.push_back(new adv_upstream<idx_ij>(nx, ny)); // x 
    advop.push_back(new adv_upstream<idx_ji>(ny, nx)); // y 
  }

  // memoray allocation: indices
  ptr_vector<rng_t> i(np), j(np);
  for (int p=0; p<np; ++p) {
    i.push_back(new rng_t(p*nx/np, (p+1)*nx/np-1));
    j.push_back(new rng_t(0, ny - 1));
  }

  // memory allocation: scalar fields at two time levels (n, n+1)
  vec_arr_t psi(ndim); 
  for (int i=0; i<adv::ntlev; ++i) 
    psi.push_back(new arr_t(advop[x].rng_psi(), advop[y].rng_psi()));

  // memory allocation: Courant number fields for two dimensions (x,y)
  vec_arr_t vel(ndim); 
  for (int i=0; i<ndim; ++i) 
    vel.push_back(new arr_t(advop[x].rng_vel(), advop[y].rng_vel()));
  vel[x] = Cx;
  vel[y] = Cy;

  // reading the initial condition (due to presence of haloes the data to be stored 
  // is not contiguous, hence looping over one dimension) dims = {t,x,y}
  for (
    pos_t ncdf_off({0,0,0}), ncdf_cnt({1,1,siz_t(ny)}); 
    ncdf_off[1] < nx; 
    ++ncdf_off[1]
  ) nv.getVar(ncdf_off, ncdf_cnt, psi[n](ncdf_off[1], rng_t(0,ny)).dataFirst());

  // integration
  for (siz_t t = 1; t <= nt; ++t) {
    for (int s = 0; s < advop[x].n_steps(); ++s) {
      // filling halos 
      for (int d = 0; d < ndim; ++d) {
        psi[n](advop[d].left_halo()) = psi[n](advop[d].rght_edge());
        psi[n](advop[d].rght_halo()) = psi[n](advop[d].left_edge());
      }

      // advecting in each dimension sharing the workload among threads
      psi[n+1] = psi[n];
#pragma omp parallel for
      for (int p = 0; p < np; ++p) // threads
        for (int d = 0; d < ndim; ++d) // dimensions
          advop[d](psi, vel[d], i[p], j[p], n, s);
    }

    // outputting
    if ((t % no) == 0) for (
      pos_t ncdf_off({siz_t(t/no),0,0}), ncdf_cnt({1,1,siz_t(ny)}); 
      ncdf_off[1] < nx; 
      ++ncdf_off[1]
    ) nv.putVar(ncdf_off, ncdf_cnt, psi[n+1](ncdf_off[1], rng_t(0,ny)).dataFirst());

    // cycling the pointers to the arrays
    blitz::cycleArrays(psi[n], psi[n+1]);
  }  
}
