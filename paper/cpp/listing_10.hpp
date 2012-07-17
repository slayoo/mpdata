template <class pi>
inline auto A(
  const arr_t &psi,
  const rng_t &i, const rng_t &j
) return_macro(
  frac(
    psi(pi(i+1, j)) - psi(pi(i,j)),
    psi(pi(i+1, j)) + psi(pi(i,j))
  ) 
) 
