//listing16
#include "listings.hpp"
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

template <class T>
void init(T &slv)
{
  // TODO: anderson_fattahi()?
  slv.state() = 0; // TODO: cone/gauss
  slv.state()(rng_t(0,2), rng_t(0,3)) = 1;  
  slv.Cx() = .25; // TODO: omega (Boost.units?) :)
  slv.Cy() = .25;
}

int main()
{
  int nx = 32, ny = 32, nt = 200;

  solver_2D<mpdata<1>, cyclic<0>, cyclic<1>> 
    slv1(nx, ny);
  solver_2D<mpdata<2>, cyclic<0>, cyclic<1>> 
    slv2(nx, ny); 
  solver_2D<mpdata<3>, cyclic<0>, cyclic<1>> 
    slv3(nx, ny); 

  Gnuplot gp;
  gp << "set term pdf size 10cm, 30cm\n";
  gp << "set output 'figure.pdf'\n";
  gp << "set nosurface\n";
  gp << "set border 4095\n";
  gp << "unset xtics\n";
  gp << "unset ytics\n";
  gp << "set ticslevel 0\n";
  //gp << "set cbrange [0:1.1]\n";
  //gp << "set zrange [0:1.1]\n";
  gp << "set multiplot layout 4,1\n";
  gp << "set pm3d\n";

  init(slv1);
  std::string fmt = gp.binfmt(slv1.state());
  gp << "set title 't=0'\n";
  gp << "splot '-' binary" << fmt << " notitle\n";
  gp.sendBinary(slv1.state().copy());

  slv1.solve(nt);
  gp << "set title 'mpdata<1> @ t=" << nt << "'\n";
  gp << "splot '-' binary" << fmt << " notitle\n";
  gp.sendBinary(slv1.state().copy());

  init(slv2); // TODO
  slv2.solve(nt);
  gp << "set title 'mpdata<2> @ t=" << nt << "'\n";
  gp << "splot '-' binary" << fmt << " notitle\n";
  gp.sendBinary(slv2.state().copy());

  init(slv3); 
  slv3.solve(nt); 
  gp << "set title 'mpdata<3> @ t=" << nt << "'\n";
  gp << "splot '-' binary" << fmt << " notitle\n";
  gp.sendBinary(slv3.state().copy());
}
//listing17
