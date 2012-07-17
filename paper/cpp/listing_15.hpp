template <class adv_t>
struct solver_2D
{
  vec_t<arr_t> psi, C;
  int n;
  rng_t i, j;
  adv_t adv;
  bcd<pi_ij, adv_t::n_halos> bcx;
  bcd<pi_ji, adv_t::n_halos> bcy;

  solver_2D(int nx, int ny) :
    n(0), i(0, nx-1), j(0, ny-1), 
    adv(i, j), bcx(i,j), bcy(j,i)
  {
    const int hlo = adv.n_halos;
    for (int l = 0; l < 2; ++l) 
      psi.push_back(new arr_t(i ^ hlo, j ^ hlo));
    C.push_back(new arr_t(i ^ h ^ hlo, j ^ hlo));
    C.push_back(new arr_t(i ^ hlo, j ^ h ^ hlo));
  }
  
  arr_t state() {return psi[n](i,j).reindex({0,0});}
  arr_t Cx() { return C[0](i ^ h, j); } // TODO: bcd!
  arr_t Cy() { return C[1](i, j ^ h); } // TODO: bcd!

  void solve(int nt) {
    for (int t = 0; t < nt; ++t) {
      for (int s = 0; s < adv.n_steps; ++s) {
        bcx.fill_halos(psi[n]);
        bcy.fill_halos(psi[n]);
        adv.op_2D(psi, n, C, i, j, s);
        n = (n + 1) % 2;
      }
    }
  }
};
