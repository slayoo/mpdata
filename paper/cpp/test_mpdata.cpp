#include "listings.hpp"
#define GNUPLOT_ENABLE_BLITZ

int main()
{
  int nx = 3, ny = 3, nt = 1;

  solver_2D<mpdata<2>, cyclic<0>, cyclic<1>> slv(nx, ny); 

  slv.state() = 0; 
  slv.state()(rng_t(1,1), rng_t(1,1)) = 1;  
  slv.Cx() = .5; 
  slv.Cy() = .5;
  slv.solve(nt);
  std::cerr << slv.state() << blitz::sum(slv.state()) << std::endl;
}
