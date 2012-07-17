struct pi_ij : blitz::RectDomain<2> 
{ 
  pi_ij(const rng_t &i, const rng_t &j) :
    RectDomain<2>(blitz::TinyVector<rng_t, 2>(i,j)) 
  {}  
};

struct pi_ji : blitz::RectDomain<2> 
{ 
  pi_ji(const rng_t &j, const rng_t &i) :
    RectDomain<2>(blitz::TinyVector<rng_t, 2>(i,j)) 
  {}  
}; 
