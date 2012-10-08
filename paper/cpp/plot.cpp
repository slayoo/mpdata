//listing16
#include "listings.hpp"
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

template <class slv_t>
void init(slv_t &slv, const int nx, const int ny)
{
  blitz::firstIndex i;
  blitz::secondIndex j;
  slv.state() = exp(
    -sqr(i-nx/2.) / (2.*pow(nx/10,2))
    -sqr(j-ny/2.) / (2.*pow(ny/10,2))
  );  
  slv.Cx() = .5; 
  slv.Cy() = .25;
}

int main()
{
  Gnuplot gp;
  gp << 
    "set term pdf size 10cm, 30cm\n"
    "set output 'figure.pdf'\n"
    "set nosurface\n"
    "set border 4095\n"
    "unset xtics\n"
    "unset ytics\n"
    "set ticslevel 0\n"
    "set cbrange [0:1.1]\n"
    "set zrange [0:1.1]\n"
    "set multiplot layout 4,1\n"
    "set pm3d\n";

  int nx = 32, ny = 32, nt = 128;
  std::string fmt;
  {
    solver_2D<mpdata<1>, cyclic<0>, cyclic<1>> 
      slv(nx, ny);
    init(slv, nx, ny);
    fmt = gp.binfmt(slv.state());
    gp << 
      "set title 't=0'\n";
      "splot '-' binary" << fmt << " notitle\n";
    gp.sendBinary(slv.state().copy());

    slv.solve(nt);
    gp << 
      "set title 'donorcell @ t=" << nt << "'\n";
      "splot '-' binary" << fmt << " notitle\n";
    gp.sendBinary(slv.state().copy());
  }
  {
    solver_2D<mpdata<2>, cyclic<0>, cyclic<1>> 
      slv(nx, ny); 
    init(slv, nx, ny); 
    slv.solve(nt);
    gp << 
      "set title 'mpdata<2> @ t=" << nt << "'\n";
      "splot '-' binary" << fmt << " notitle\n";
    gp.sendBinary(slv.state().copy());
  }
  {
    solver_2D<mpdata<3>, cyclic<0>, cyclic<1>> 
      slv(nx, ny); 
    init(slv, nx, ny); 
    slv.solve(nt); 
    gp << 
      "set title 'mpdata<3> @ t=" << nt << "'\n";
      "splot '-' binary" << fmt << " notitle\n";
    gp.sendBinary(slv.state().copy());
  }
}
//listing17
