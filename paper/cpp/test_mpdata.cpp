#include "listings.hpp"
int main()
{
  int nx = 3, ny = 3, nt = 1;

  {
    solver_2D<mpdata<1>, cyclic<0>, cyclic<1>> slv(nx, ny); 

    slv.state() = 0;  
    slv.state()(1,1) = 1;  
    slv.Cx() = .2; 
    slv.Cy() = .2; 

    std::cerr << slv.state() << std::endl;
    slv.solve(nt);
    std::cerr << slv.state() << std::endl;
  }
  {
    solver_2D<mpdata<2>, cyclic<0>, cyclic<1>> slv(nx, ny); 
    slv.state() = 0;  
    slv.state()(1,1) = 1;  
    slv.Cx() = .2; 
    slv.Cy() = .2; 
    slv.solve(nt);
    std::cerr << slv.state() << std::endl;
  }
  {
    solver_2D<mpdata<3>, cyclic<0>, cyclic<1>> slv(nx, ny); 
    slv.state() = 0;  
    slv.state()(1,1) = 1;  
    slv.Cx() = .2; 
    slv.Cy() = .2; 
    slv.solve(nt);
    std::cerr << slv.state() << std::endl;
  }
}
