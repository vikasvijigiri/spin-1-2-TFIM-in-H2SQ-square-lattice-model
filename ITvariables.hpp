#ifndef _ITVARIABLES_HPP_DEFINED_
#define _ITVARIABLES_HPP_DEFINED_
#include <random>
#include <ctime>
#pragma once

extern	int lx, ly, Ns, Nl, Np, M;
extern	double J0,J1,J2,J3,J3_temp,h,Kx,T;
extern	int nbins, isteps, iter;
extern  std::string coupling;
extern  std::string versus;
extern bool bool_J3;

extern  int myrank;
extern  int numprocs;
extern  int exc_index;

class ITvariables {
	public:
      	void declare_variables(int );
      	void set_temperatures(int );

};

class ran {
  private:
    std::mt19937 mt;
    std::uniform_real_distribution < double > dist;

  public:
    ran(double lower = 0.0, double upper = 1.0): mt(std::time(nullptr)), dist(lower, upper) {}
    double operator()() {
    return dist(mt);
  }
};
#endif
