#include "listings.hpp"
int main()
{
  int nx = 3, ny = 3, nt = 1;

  solver_2D<mpdata<2>, cyclic<0>, cyclic<1>> slv(nx, ny); 

  slv.state() = 0; 
  slv.state()(1,1) = 1;  
  slv.Cx() = .5; 
  slv.Cy() = .1;

  //if (blitz::sum(slv.state()) != 1) throw;
  //else std::cerr << "sum OK" << std::endl;

  std::cerr << slv.state() << std::endl;
  slv.solve(nt);
  std::cerr << slv.state() << std::endl;

/*
  std::cerr << "checking values..." << std::endl;
  if (
    slv.state()(0,0) != .5 ||
    slv.state()(1,0) != .5 ||
    slv.state()(0,1) != 0 ||
    slv.state()(1,1) != 0
  ) throw;
  std::cerr << "OK" << std::endl;
*/
}
