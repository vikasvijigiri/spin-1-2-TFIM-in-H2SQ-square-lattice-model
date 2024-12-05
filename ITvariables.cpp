#include "ITlattice.hpp"
#include "ITvariables.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>

int lx, ly, Ns, Nl, Np, M;
double J0, J1, J2, J3, J3_temp, Kx, h, T;
int nbins, isteps, iter;
std::string coupling;
std::string versus;
int myrank, numprocs, exc_index;
bool bool_J3;

void
ITvariables::declare_variables(int myrank)
{
	// Declares initial mandatory parameters.
	std::ifstream vars("/data/docs/vikas/IT_QMC_H2SQ/input_params.in");

	vars >> lx;
	vars >> J0;
	vars >> J1;
	vars >> J2;
	vars >> J3;
	vars >> M;
	vars >> Nl;
	vars >> nbins;
	vars >> iter;
	vars >> isteps;
	ly = lx;
	Np = lx * ly;
	Ns = 2 *lx * ly;	//  Simple square lattice;
	if (J3 == 0)
	{
		bool_J3 = false;
	}
	else {
		bool_J3 = true;
		J3_temp=J3;
		//J3=0;
	}
	coupling = "AFM";
	versus = "T";
}

void
ITvariables::set_temperatures(int myrank)
{
	if (versus == "T")
	{
		// Only temperature is varying;
		T = (myrank+1) *0.02f;
		h = 0.f;
		Kx = 0.f;
	}
	else if (versus == "TF")
	{
		// Both temperature and field are varying according to theta = arcTan(beta*Kx);
		T = 0.05f + (myrank) *0.05f;
		double Tn = 6.0f;
		double theta = (atan(1.f) *4.f / Tn);
		h = T* tan(theta);
		Kx = -T *((log(tanh(h / (T *M)))) / 2.f);
	}
	else if (versus == "F")
	{
		// Onlyfield is varying;
		h = 0.05f + (myrank) *0.1f;
		T = 0.02f;
		Kx = -T *((log(tanh(h / (T *M)))) / 2.f);
	}
}
