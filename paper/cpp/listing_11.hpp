template <class pi>
inline auto B(
  const arr_t &psi, 
  const rng_t &i, const rng_t &j
) return_macro(
 .5 * frac(
    psi(pi(i+1, j+1)) + psi(pi(i, j+1)) -
    psi(pi(i+1, j-1)) - psi(pi(i, j-1)),
    psi(pi(i+1, j+1)) + psi(pi(i, j+1)) +
    psi(pi(i+1, j-1)) + psi(pi(i, j-1))
  )
)
