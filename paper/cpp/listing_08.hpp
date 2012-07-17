inline void donorcell_2D(
  const vec_t<arr_t> &psi, const int n,
  const vec_t<arr_t> &C, 
  const rng_t &i, const rng_t &j
) { 
  psi[n+1](i,j) = psi[n](i,j)
    - donorcell_1D<pi_ij>(psi[n], C[0], i, j)
    - donorcell_1D<pi_ji>(psi[n], C[1], j, i); 
}
