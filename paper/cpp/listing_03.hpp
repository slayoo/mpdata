struct hlf_t {} h;

inline rng_t operator+(
  const rng_t &i, const hlf_t &
) { 
  return i + 1; 
} 

inline rng_t operator-(
  const rng_t &i, const hlf_t &
) { 
  return i; 
} 
