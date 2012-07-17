template <int n_iters>
struct mpdata 
{
  // member fields
  vec_t<arr_t> tmp0, tmp1;
  static const int n_steps = n_iters;
  static const int n_halos = n_iters;

  // ctor
  mpdata(const rng_t &i, const rng_t &j) 
  {
    const int hlo = n_halos;
    tmp0.push_back(new arr_t(i ^ h ^ hlo, j ^ hlo));
    tmp0.push_back(new arr_t(i ^ hlo, j ^ h ^ hlo));
    if (n_iters < 2) return;
    tmp1.push_back(new arr_t(i ^ h ^ hlo, j ^ hlo));
    tmp1.push_back(new arr_t(i ^ hlo, j ^ h ^ hlo));
  }

  // the 2D advection operator method
  inline void op_2D( 
    const vec_t<arr_t> &psi, const int n,
    const vec_t<arr_t> &C, 
    const rng_t &i, const rng_t &j,
    const int step
  ) {
    // chosing input/output for antidiff. velocity
    const vec_t<arr_t> 
      &C_unco = (step == 1) ? C : 
        (step % 2) ? tmp0 : tmp1,
      &C_corr = (step  % 2) ? tmp1 : tmp0;
    
    // calculating the antidiffusive velocities
    if (step > 0) {
      const int hlo = n_steps - 1 - step;
      enum {x,y};
      C_corr[x](i ^ h ^ hlo, j ^ hlo) =
        antidiff_2D<pi_ij>(
          psi[n], i ^ h ^ hlo, j ^ hlo, C_unco, x
        );
      C_corr[y](i ^ hlo, j ^ h ^ hlo) =
        antidiff_2D<pi_ji>(
          psi[n], j ^ h ^ hlo, i ^ hlo, C_unco, y
        );
    }

    // performing a donor-cell step with C or C_corr
    donorcell_2D(psi, n, step==0 ? C : C_corr, i, j);
  }
};
