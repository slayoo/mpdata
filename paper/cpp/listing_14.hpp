template <class pi, int hlo>
struct bcd
{
  // (could-be-private) member fields
  pi left_halo, rght_edge, rght_halo, left_edge;

  // ctor with member-field initialisation 
  bcd(const rng_t &i, const rng_t &j) :
    //     ( i.first  ...   i.last         ), (  j ) 
    left_halo(
      rng_t(i.first()-hlo,  i.first()-1    ), j^hlo),
    rght_edge(
      rng_t(i.last()-hlo+1, i.last()       ), j^hlo),
    rght_halo(
      rng_t(i.last()+1,     i.last()+hlo   ), j^hlo),
    left_edge(
      rng_t(i.first(),      i.first()+hlo-1), j^hlo)
  {} 

  // public method
  void fill_halos(const arr_t &psi)
  {
    psi(l_halo) = psi(r_edge);     
    psi(r_halo) = psi(l_edge);     
  }
};
