template <class nom_t, class den_t>
static inline auto frac(
  const nom_t &nom, const den_t &den
) return_macro(
  where(den > 0, nom / den, 0)
) 
