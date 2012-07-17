template <class pi>  
inline auto donorcell_1D( 
  const arr_t &psi, 
  const arr_t &C, 
  const rng_t &i, 
  const rng_t &j
) return_macro(
  F(psi(pi(i,  j)), psi(pi(i+1,j)), C(pi(i+h,j))) -
  F(psi(pi(i-1,j)), psi(pi(i,  j)), C(pi(i-h,j)))
)
