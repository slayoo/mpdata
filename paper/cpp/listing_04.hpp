template <class n_t>
inline rng_t operator^(
  const rng_t &r, const n_t &n) 
{ 
  return rng_t(
    (r - n).first(), 
    (r + n).last()
  ); 
} 
