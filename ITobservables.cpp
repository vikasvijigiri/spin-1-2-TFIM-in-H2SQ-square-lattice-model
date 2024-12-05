#include "ITlattice.hpp"
#include "ITvariables.hpp"
#include "ITupdates.hpp"
#include "ITobservables.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>

void
ITobservables::observables(ITlattice ***lattice, ITupdates &upd)
{

	double k_x = 0.f, k_y = 3.14159265359f, Icl,x,y,p;
	double ms1 = 0.f, ms2 = 0.f, mc1 = 0.f, mc2 = 0.f, cs = 1.f;
	double Enrg = 0.f;
	int zc = 0;
	for (int xc = 0; xc < M; ++xc)
	{
		for (int yc = 0; yc < Ns; ++yc)
		{
			//!-----------------------Energy and specific heat---------------------------!
			for (int zcR = 0; zcR < Nl; ++zcR)
			{
			  Enrg+=upd.enrg_diff_local(lattice, xc, yc, zcR, 1)/2;
			}

			//!------------------ stripe order P, rho (X_P, X_rho)-----------------------!
			double w = cs *((k_x) *psite[yc][0] + (k_y) *psite[yc][1]);
			ms1 += lattice[xc][yc][zc].S() *sin(w);
			mc1 += lattice[xc][yc][zc].S() *cos(w);
			w = cs *((k_y) *psite[yc][0] + (k_x) *psite[yc][1]);
			ms2 += lattice[xc][yc][zc].S() *sin(w);
			mc2 += lattice[xc][yc][zc].S() *cos(w);
			//!------------------------------------------------------!
			if (yc < Np)
			{
				int v1 = lattice[xc][yc][zc].S() - lattice[xc][upd.PBC(yc, lx, Np)][zc].S() + lattice[xc][yc + Np][zc].S() - lattice[xc][upd.PBC(yc, Np + 1, lx)][zc].S();
				int v2 = lattice[xc][yc][zc].S() - lattice[xc][upd.PBC(yc, lx, Np)][zc].S() - lattice[xc][yc + Np][zc].S() + lattice[xc][upd.PBC(yc, Np + 1, lx)][zc].S();
				if (abs(v1) == 4)
				{
					if (yc < lx *ly)
					{
						Icl += 1;
					}
				}
				else if (abs(v2) == 4)
				{
					if (yc < lx *ly)
					{
						Icl += 1;
					}
				}
				else if (abs(v1) != 4 or abs(v2) != 4)
				{
					if (yc < lx *ly)
					{
						Icl -= 1 / 3;
					}
				}
			}
		}
	}
	x=pow(pow(pow(ms1, 2) + pow(mc1, 2), 2) + pow(pow(ms2, 2) + pow(mc2, 2), 2), 0.5)/(pow(Ns * 1 * M, 2));
	op += x;
	y=Icl/(lx * ly * 1 * M);
	Ic += y;
	E += Enrg/(Ns * Nl * M);
	/*
	   	// Correlation function in imaginary-time, to detect deconfined phase.
	    if (J2 == 0 and M >=32) {
	        for (int xc = 10; xc < 32; ++xc) 
	        {
	            for (int yc = 0; yc<lx*ly; ++yc) 
	            { 
	            	//for (int zc = 0; zc < ly; ++zc) 
			//{
	                     double P1 = lattice[0][yc].S();//-(1/4)*(lattice[0][yc].S() - lattice[0][PBC(yc,1,lx)].S() + lattice[0][yc+lx].S() - lattice[0][yc+lx][PBC(zc,1,lx)].S());
	                     double P2 = lattice[xc][yc].S();// -(1/4)*(lattice[xc][yc].S() - lattice[xc][PBC(yc,1,lx)].S() + lattice[xc][yc+lx].S() - lattice[xc][yc+lx][PBC(zc,1,lx)].S());
	                     c[xc - 10] += P1 *P2 / (lx*ly);
	               	//}
	            }
	        }
	             	// std::cout << "vikas" <<std::endl;     
	    *enrg = *enrg / (2 *lx *ly *Nl);
	    *enrg_sq = pow(*enrg, 2);
	    *m1 = abs(*m1 / (M *2 *lx *ly *Nl));
	    *m1_sq = pow(*m1, 2);
	    */
	op_sq += pow(x, 2);
	op_four += pow(x, 4);
	Ic_sq += (y) *(y);
	E_sq += pow(Enrg, 2);	

}
void
ITobservables::binning_data(int bins, int iters)
{

	op /= iters;
	op_sq /= iters;
	op_four /= iters;
	Ic /= iters;
	Ic_sq /= iters;
	E /= iters;
	E_sq /= iters;


	data1[0] += op;
	data1[1] += op_sq;
	data1[2] += op_four;
	data1[3] += Ic;
	data1[4] += Ic_sq;
	data1[5] += E;
	data1[6] += E_sq;	

	data2[0] += pow(op, 2);
	data2[1] += pow(op_sq, 2);
	data2[2] += pow(op_four, 2);
	data2[3] += pow(Ic, 2);
	data2[4] += pow(Ic_sq, 2);
	data2[5] += pow(E, 2);
	data2[6] += pow(E_sq, 2);	

	for (int i = 0; i < 7; ++i)
	{
		wdata1[i] = data1[i] / bins;
		wdata2[i] = data2[i] / bins;
		wdata2[i] = sqrt(abs(wdata2[i] - pow(wdata1[i], 2)) / bins);
	}

	op = 0.;
	op_sq = 0.;
	op_four = 0.;
	Ic = 0.;
	Ic_sq = 0.;
	E=0.;
	E_sq=0.;	
}

void
ITobservables::Initiate_observables()
{
	data1 = new double[7];
	data2 = new double[7];
	wdata1 = new double[7];
	wdata2 = new double[7];
	op=0; op_sq=0; op_four=0; Ic=0; Ic_sq=0; E=0; E_sq=0;
}

void ITobservables::pbc()
{
	int k, x1, y1, d1, d2;
	psite = new int *[Ns];
	for (int i = 0; i < Ns; i++) psite[i] = new int[2];
	for (k = 1; k <= 2 *lx * ly; ++k)
	{
		// Block 1(Bottom left)
		if (k <= lx *ly)
		{
			if (k % lx <= lx / 2 && k <= (int((k - 1) / lx) + 1) *lx && k <= (Np / 2) + lx)
			{
				if (k % lx != 0)
				{
					x1 = (k % lx) - 1;
					y1 = int((k - 1) / lx);
					d1 = x1 + y1;
					d2 = y1 - x1;
				}
				else
				{
					x1 = (-k % lx) - 1;
					y1 = int((k - 1) / lx);
					d1 = x1 + y1;
					d2 = y1 - x1;
				}
			}
			// Block 2(Botton right)
			else if (k % lx > lx / 2 && k <= (int((k - 1) / lx) + 1) *lx && k <= (Np / 2) + lx)
			{
				x1 = (k % lx) - 1 - lx;
				y1 = int((k - 1) / lx);
				d1 = x1 + y1;
				d2 = y1 - x1;
			}
			// Block 3(top left)
			else if (k % lx <= lx / 2 && k <= (int((k - 1) / lx) + 1) *lx && k >= (Np / 2) + lx)
			{
				x1 = (k % lx) - 1;
				y1 = int((k - 1) / lx) - lx;
				d1 = x1 + y1;
				d2 = y1 - x1;
			}
			// Block 4(top right)
			else if (k % lx > lx / 2 && k <= (int((k - 1) / lx) + 1) *lx && k >= (Np / 2) + lx)
			{
				x1 = (k % lx) - 1 - lx;
				y1 = int((k - 1) / lx) - lx;
				d1 = x1 + y1;
				d2 = y1 - x1;
			}
			psite[k - 1][0] = d1;
			psite[k - 1][1] = d2;
		}
		else
		{
			if ((k - 1) % lx <= lx / 2 && k <= (int((k - 1) / lx) + 1) *lx && k <= 2 *Np - Np / 2)
			{
				x1 = k % lx - 1 + int(k / lx) - lx;
				y1 = int(((k - Np) - 1) / lx) + 1 - (k - 1) % lx;
			}
			else if ((k - 1) % lx > lx / 2 && k <= (int((k - 1) / lx) + 1) *lx && k <= 2 *Np - Np / 2)
			{
				x1 = (k - 1) % lx - lx + int((k - Np - 1) / lx);
				y1 = lx - (k - 1) % lx + int((k - Np - 1) / lx) + 1;
			}
			else if ((k - 1) % lx <= lx / 2 && k <= (int((k - 1) / lx) + 1) *lx && k >= 2 *Np - Np / 2)
			{
				x1 = (k - Np) % lx - 1 + int((k - Np) / lx) - lx;
				y1 = int(((k - 2 *Np)) / lx) - (k - Np - 1) % lx;
			}
			else if ((k - 1) % lx > lx / 2 && k <= (int((k - 1) / lx) + 1) *lx && k >= 2 *Np - Np / 2)
			{
				x1 = (k - 1) % lx - lx - int((2 *Np - k + 1) / lx) - 1;
				y1 = -(k - 1) % lx + lx - int((2 *Np - k + 1) / lx);
			}

			psite[k - 1][0] = x1;
			psite[k - 1][1] = y1;
		}
	}	// for loop
}	// end func pbc
