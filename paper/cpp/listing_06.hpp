template <class T1, class T2, class T3> 
inline auto F(
  const T1 &psi_l, 
  const T2 &psi_r, 
  const T3 &C
) return_macro(
  .5 * (C + abs(C)) * psi_l + 
  .5 * (C - abs(C)) * psi_r
)
