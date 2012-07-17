template <class pi>
inline auto antidiff_2D(
  const arr_t &psi, 
  const rng_t &i, const rng_t &j,
  const vec_t<arr_t> &C, const int d
) return_macro(
  abs(C[d](pi(i+h, j)))
  * A<pi>(psi, i, j)
  * (1 - abs(C[d](pi(i+h, j))))
  - C[d](pi(i+h, j)) 
  * B<pi>(psi, i, j)
  * .25 * (
    C[d+1](pi(i+1, j+h)) + C[d+1](pi(i, j+h)) +
    C[d+1](pi(i+1, j-h)) + C[d+1](pi(i, j-h)) 
  ) 
) 
