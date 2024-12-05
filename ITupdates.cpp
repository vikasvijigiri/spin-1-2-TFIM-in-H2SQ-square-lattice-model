#include "ITlattice.hpp"
#include "ITvariables.hpp"
#include "ITupdates.hpp"
#include "ITobservables.hpp"
#include <fstream>
#include <iostream>
#include <deque>
#include <iomanip>
#include <cstdlib>

static ran rann;
using namespace std;

void ITupdates::make_bonds(ITlattice ***lattice, ITlattice *plaquette, deque<int> *cluster, deque<int> *plaq_cluster, int xc, int yc, int zc, double temperature)
{
	//double beta = 1 / temperature;
	ITlattice central = lattice[xc][yc][zc];
	ITlattice plaq = plaquette[yc % Np];
	int ycR = yc;
	int n1, n2, n3, n4, n5, n6;

	for (int i = 0; i < 2; ++i)
	{
		yc = ycR;
		while (1)
		{
			int nbr = 0;
			if (yc < Np)
			{
				n1 = PBXR(yc + Np, 1, lx);    
				n2 = PBYU(yc, lx, Np); 
				n3 = yc + Np;
				n4 = PBYD(n3, -lx, Np);
				n5 = PBYD(yc, -lx, Np);
				n6 = PBXR(n4, 1, lx);
			}
			else
			{
				n1 = yc % Np; 
				n2 = PBXR(yc, 1, lx); 
				n3 = PBYU(n1, lx, Np);
				n4 = PBXL(n3, -1, lx);
				n5 = PBXL(yc, -1, lx);
				n6 = PBYD(n4, -lx, Np);
			}	

			int G1 = lattice[xc][yc][zc].S() *lattice[xc][n1][zc].S() *lattice[xc][n2][zc].S() *lattice[xc][n3][zc].S();
			int G2 = lattice[xc][yc][zc].S() *lattice[xc][n4][zc].S() *lattice[xc][n5][zc].S() *lattice[xc][n6][zc].S();
			if (rann() < central.ice_rule(lattice[xc][n2][zc]) and G1 > 0)
			{
			  if (!lattice[xc][n2][zc].bonded() and !plaquette[yc % Np].bonded()) {
					lattice[xc][n2][zc].set_cl(central.cl());
					plaquette[yc % Np].set_cl(plaq.cl());
					plaq = plaquette[yc % Np];
					nbr = 1;
					if (n2 != ycR) cluster->push_back(n2);
					plaq_cluster->push_back(yc % Np);
					yc = n2;
			  }		
			}
			else if (rann() < central.ice_rule(lattice[xc][n5][zc]) and G2 > 0)
			{
			  if (!lattice[xc][n5][zc].bonded() and !plaquette[n5 % Np].bonded()){
					lattice[xc][n5][zc].set_cl(central.cl());
					plaquette[n5 % Np].set_cl(plaq.cl());
					plaq = plaquette[n5 % Np];
					nbr = 1;
					if (n5 != ycR) cluster->push_back(n5);
					plaq_cluster->push_back(n5 % Np);
					yc = n5;
			  }		
			}
			else if (rann() < central.ice_rule(lattice[xc][n1][zc]) and G1 > 0)
			{
			  if (!lattice[xc][n1][zc].bonded() and !plaquette[yc % Np].bonded()){
					lattice[xc][n1][zc].set_cl(central.cl());
					plaquette[yc % Np].set_cl(plaq.cl());
					plaq = plaquette[yc % Np];
					nbr = 1;
					if (n1 != ycR) cluster->push_back(n1);
					plaq_cluster->push_back(yc % Np);
					yc = n1;
			  }		
			}
			else if (rann() < central.ice_rule(lattice[xc][n3][zc]) and G1 > 0)
			{
			  if (!lattice[xc][n3][zc].bonded() and !plaquette[yc % Np].bonded() ){
					lattice[xc][n3][zc].set_cl(central.cl());
					plaquette[yc % Np].set_cl(plaq.cl());
					plaq = plaquette[yc % Np];
					nbr = 1;
					if (n3 != ycR) cluster->push_back(n3);
					plaq_cluster->push_back(yc % Np);
					yc = n3;
			  }		
			}
			else if (rann() < central.ice_rule(lattice[xc][n4][zc]) and G2 > 0)
			{
			 if (yc < Np){
			  if (!lattice[xc][n4][zc].bonded() and !plaquette[n5 % Np].bonded()){
			  		plaquette[n5 % Np].set_cl(plaq.cl()); plaq_cluster->push_back(n5 % Np);	plaq = plaquette[n5 % Np];
					lattice[xc][n4][zc].set_cl(central.cl());			  
					nbr = 1;  
					if (n4 != ycR) cluster->push_back(n4);
					yc = n4;													
			  }
			 }
			 else  {       
			  if (!lattice[xc][n4][zc].bonded() and !plaquette[n6 % Np].bonded()){			  
			  		plaquette[n6 % Np].set_cl(plaq.cl()); plaq_cluster->push_back(n6 % Np);	plaq = plaquette[n6 % Np];
					lattice[xc][n4][zc].set_cl(central.cl());			  
					nbr = 1;
					if (n4 != ycR) cluster->push_back(n4);
					yc = n4;
			  }
			 }	
			}
			else if (rann() < central.ice_rule(lattice[xc][n6][zc]) and G2 > 0)
			{
			 if (yc < Np){
			  if (!lattice[xc][n6][zc].bonded() and !plaquette[n5 % Np].bonded()){
			 		plaquette[n5 % Np].set_cl(plaq.cl()); plaq_cluster->push_back(n5 % Np);	plaq = plaquette[n5 % Np];
					lattice[xc][n6][zc].set_cl(central.cl());
					nbr = 1;
					if (n6 != ycR) cluster->push_back(n6);
					yc = n6;
			  }
			 } 		
			 else {
			  if (!lattice[xc][n6][zc].bonded() and !plaquette[n6 % Np].bonded()){
			 		lattice[xc][n6][zc].set_cl(central.cl());
					plaquette[n6 % Np].set_cl(plaq.cl()); plaq_cluster->push_back(n6 % Np);	plaq = plaquette[n6 % Np];
					nbr = 1;
					if (n6 != ycR) cluster->push_back(n6);
					yc = n6;
			  }
			 } 		
			}
			central = lattice[xc][yc][zc];
			//**************************************************************************************************************************************

			if (ycR == yc or nbr == 0) break;
		}
	}	
}









void ITupdates::flip_cluster(ITlattice ***lattice, deque<int> cluster, bool decision, int xc, int zc, double temperature)
{
	//double beta=1/temperature;
	if (!decision)
	{
		for (std::deque<int>::iterator it = cluster.begin(); it != cluster.end(); ++it)
		{
			//for (int zc=zcR; zc < Nl; ++zc){		//    	ITlattice central = lattice[xc][*it][zc];
			lattice[xc][*it][zc].flip();
			lattice[xc][*it][zc].reset_bond();		//	Here bond is like site.
			//    if (rann() > central.bond_prob(lattice[xc][*it][zc], beta, J3 / Ml))break;
			//}
			//for (int zc=zcR; zc >= 0; --zc){		//    	ITlattice central = lattice[xc][*it][zc];
			//    lattice[xc][*it][zc].flip();
			//    lattice[xc][*it][zc].reset_bond();	//	Here bond is like site.
			//    if (rann() > central.bond_prob(lattice[xc][*it][zc], beta, J3 / Ml))break;
			//}
		}
	}
}

void ITupdates::merge_clusters_intl(ITlattice ***lattice, int cluster_id1, int cluster_id2)
{
	int merged_cluster_id = std::min(cluster_id1, cluster_id2);
	int max_cluster_id = std::max(cluster_id1, cluster_id2);
	for (int xc = 0; xc < M; xc++)
	{
		for (int yc = 0; yc < Ns; yc++)
		{
			for (int zc = 0; zc < Nl; zc++)
			{
				if (lattice[xc][yc][zc].cl() == max_cluster_id) lattice[xc][yc][zc].set_cl(merged_cluster_id);
			}
		}
	}
}

void ITupdates::make_bonds_intl(ITlattice ***lattice, int xc, int yc, int zc, double temperature)
{
	double beta = 1 / temperature;
	ITlattice central = lattice[xc][yc][zc];

	if (rann() < central.bond_prob(lattice[xc][yc][PBC(zc, 1, Nl)], beta, J3 / Ml))
	{
		if (!lattice[xc][yc][PBC(zc, 1, Nl)].bonded())
		{
			lattice[xc][yc][PBC(zc, 1, Nl)].set_cl(central.cl());
		}
		else
		{
			merge_clusters_intl(lattice, lattice[xc][yc][PBC(zc, 1, Nl)].cl(), central.cl());
		}
	}

	if (rann() < central.bond_prob(lattice[xc][yc][PBC(zc, -1, Nl)], beta, J3 / Ml))
	{
		if (!lattice[xc][yc][PBC(zc, -1, Nl)].bonded())
		{
			lattice[xc][yc][PBC(zc, -1, Nl)].set_cl(central.cl());
		}
		else
		{
			merge_clusters_intl(lattice, lattice[xc][yc][PBC(zc, -1, Nl)].cl(), central.cl());
		}
	}

	if (rann() < central.bond_prob(lattice[PBC(xc, 1, M)][yc][zc], beta, Kx))
	{
		if (!lattice[PBC(xc, 1, M)][yc][zc].bonded())
		{
			lattice[PBC(xc, 1, M)][yc][zc].set_cl(central.cl());
		}
		else
		{
			merge_clusters_intl(lattice, lattice[PBC(xc, 1, M)][yc][zc].cl(), central.cl());
		}
	}

	if (rann() < central.bond_prob(lattice[PBC(xc, -1, M)][yc][zc], beta, Kx))
	{
		if (!lattice[PBC(xc, -1, M)][yc][zc].bonded())
		{
			lattice[PBC(xc, -1, M)][yc][zc].set_cl(central.cl());
		}
		else
		{
			merge_clusters_intl(lattice, lattice[PBC(xc, -1, M)][yc][zc].cl(), central.cl());
		}
	}
}

int ITupdates::make_cluster_intl(ITlattice ***lattice, int yc, double temperature)
{
	int cluster_id = 0;
	for (int xc = 0; xc < M; xc++)
	{
		for (int zc = 0; zc < Nl; zc++)
		{
			if (!lattice[xc][yc][zc].bonded())
			{
				lattice[xc][yc][zc].set_cl(cluster_id);
			}
			make_bonds_intl(lattice, xc, yc, zc, temperature);

			if (lattice[xc][yc][zc].cl() == cluster_id)
			{
				cluster_id++;
			}
		}
	}
	return cluster_id;
}

void ITupdates::flip_cluster_intl(ITlattice ***lattice, int yc, bool *decision_list)
{
	for (int xc = 0; xc < M; xc++)
	{
		for (int zc = 0; zc < Nl; zc++)
		{
			if (decision_list[lattice[xc][yc][zc].cl()]) lattice[xc][yc][zc].flip();
			lattice[xc][yc][zc].reset_bond();
		}
	}
}

double *ITupdates::cluster_field_intl(ITlattice ***lattice, int no_clusters, int yc)
{
	double *cluster_field = new double[no_clusters];
	for (int xc = 0; xc < M; xc++)
	{
		for (int zc = 0; zc < Nl; zc++)
		{
			cluster_field[lattice[xc][yc][zc].cl()] += enrg_diff_local(lattice, xc, yc, zc, 0, 0);
		}
	}
	return cluster_field;
}



int ITupdates::make_cluster_loop(ITlattice ***lattice, deque<int> *cluster, double temperature)
{
	int cluster_id = 0;
	for (int xc = 0; xc < M; xc++)
	{
		for (int zc = 0; zc < Nl; zc++)
		{
		   	for (std::deque<int>::iterator it = cluster.begin(); it != cluster.end(); ++it)
			{
				int yc=*it;
				if (!lattice[xc][yc][zc].bonded())
				{
					lattice[xc][yc][zc].set_cl(cluster_id);
				}
				make_bonds_loop(lattice, xc, yc, zc, temperature);

				if (lattice[xc][yc][zc].cl() == cluster_id)
				{
					cluster_id++;
				}
			}
		}
	}
	return cluster_id;
}

void ITupdates::flip_cluster_loop(ITlattice ***lattice, deque<int> *cluster, bool *decision_list)
{
	for (int xc = 0; xc < M; xc++)
	{
		for (int zc = 0; zc < Nl; zc++)
		{
		   	for (std::deque<int>::iterator it = cluster.begin(); it != cluster.end(); ++it)
			{
				int yc=*it;		    
				if (decision_list[lattice[xc][yc][zc].cl()]) lattice[xc][yc][zc].flip();
				lattice[xc][yc][zc].reset_bond();
			}	
		}
	}
}

double *ITupdates::cluster_field_loop(ITlattice ***lattice, int no_clusters, deque<int> *cluster)
{
	double *cluster_field = new double[no_clusters];
	for (int xc = 0; xc < M; xc++)
	{
		for (int zc = 0; zc < Nl; zc++)
		{
		   	for (std::deque<int>::iterator it = cluster.begin(); it != cluster.end(); ++it)
			{
				int yc=*it;		
				cluster_field[lattice[xc][yc][zc].cl()] += enrg_diff_local(lattice, xc, yc, zc, 0, 0);
			}	
		}
	}
	return cluster_field;
}








void ITupdates::swendsen_flip(ITlattice ***lattice, ITlattice *plaquette, double temperature)
{

	// ITlattice Update both along the interlayers and trotter axis. 
	int no_clusters;
	for (int ycR = 0; ycR < Ns; ++ycR)
	{	int yc = rann() *Ns; 
		no_clusters = make_cluster_intl(lattice, yc, temperature);
		bool *flip_decision = new bool[no_clusters];
		double *c_dE = new double[no_clusters];
		c_dE = cluster_field_intl(lattice, no_clusters, yc);
		for (int i = 0; i < no_clusters; i++) flip_decision[i] = (rann() <= exp(c_dE[i] / temperature) ); //
		flip_cluster_intl(lattice, yc, flip_decision);
		delete[] flip_decision;
		delete[] c_dE;
	}
	
	// Loop Update in 2D along with the single flip metropolis in all directions.  
	double dE,r; int yc;	
	deque<int> cluster, plaq_cluster;
	
	for (int ycR = 0; ycR < Ns; ++ycR)
	{
		int yc = rann()*Ns;
		cluster.push_back(yc);
		make_bonds(lattice, plaquette, &cluster, &plaq_cluster, xc, yc, zc, temperature);		// filled cluster
		
		no_clusters = make_cluster_intl(lattice, cluster, temperature);
		bool *flip_decision = new bool[no_clusters];
		double *c_dE = new double[no_clusters];
		c_dE = cluster_field_loop(lattice, no_clusters, cluster);
		for (int i = 0; i < no_clusters; i++) flip_decision[i] = (rann() <= exp(c_dE[i] / temperature) ); //
		flip_cluster_loop(lattice, cluster, flip_decision);
		delete[] flip_decision;
		delete[] c_dE;

		plaq_cluster = dip_enrg_diff_cluster(lattice, plaquette, cluster, plaq_cluster, &dE, xc, zc);
		flip_cluster(lattice, cluster, (r <= exp(dE / temperature)), xc, zc, temperature);	// perform operations and empty the cluster and return back.
		cluster.clear();
	}

	
}
////////////////////////////////////////////////////////////////////////////////////////////////////


double ITupdates::enrg(ITlattice ***lattice, int J3_prefac, int K_prefac)
{
	E=0;
	double px1, py1, px[2], py[2];
	int n[2];
	for (int xc = 0; xc < M; ++xc)
	{
		for (int zc = 0; zc < Nl; ++zc)
		{
			for (int yc = 0; yc < Np; ++yc)
			{
				n[0] = PBYU(yc, lx, Np);
				n[1] = PBXL(yc, -1, lx);
				for (int i = 0; i < 2; ++i)
					dip(lattice, xc, n[i], zc, &px[i], &py[i]);
				E += -(J0 / Ml) *(lattice[xc][yc][zc].S() *lattice[xc][PBYU(yc, lx, Np)][zc].S() *lattice[xc][yc + Np][zc].S() *lattice[xc][PBXR(yc + Np, 1, lx)][zc].S()) +
				(J1 / Ml) *(lattice[xc][yc][zc].S() *lattice[xc][PBYU(yc, lx, Np)][zc].S() + lattice[xc][yc + Np][zc].S() *lattice[xc][PBXR(yc + Np, 1, lx)][zc].S()) +
				(J2 / Ml) *(px1*(px[0]+px[1])+py1*(py[0]+py[1])) +
				J3_prefac *(J3 / Ml) *(lattice[xc][yc][PBC(zc, 1, Nl)].S() * lattice[xc][yc][zc].S()) -
				K_prefac * Kx *(lattice[PBC(xc, 1, M)][yc][zc].S() * lattice[xc][yc][zc].S())
				;
			}
		}
	}
	return E;

}

double ITupdates::enrg_diff_local(ITlattice ***lattice, int xc, int yc, int zc, int J3_prefac, int K_prefac)
{

	double dE, p;
	/*
	dE = 2 *lattice[xc][yc][zc].S() *(-(J0 / Ml) *(lattice[xc][PBYU(yc, lx, Np)][zc].S() + lattice[xc][PBYD(yc, -lx, Np)][zc].S()+
						       lattice[xc][PBXR(yc, 1, lx)][zc].S()  + lattice[xc][PBXL(yc, -1, lx)][zc].S()) 
		- K_prefac * Kx *(lattice[PBC(xc, 1, M)][yc][zc].S() + lattice[PBC(xc, -1, M)][yc][zc].S())
	);
	*/
	
	pl(lattice, xc, yc, zc, &p);

	if (yc < Np)
	{
		dE = 2 *lattice[xc][yc][zc].S() *(-(J0 / Ml) *(lattice[xc][PBYU(yc, lx, Np)][zc].S() *lattice[xc][yc + Np][zc].S() *lattice[xc][PBXR(yc + Np, 1, lx)][zc].S() +
				lattice[xc][PBYD(yc, -lx, Np)][zc].S() *lattice[xc][PBYD(yc + Np, -lx, Np)][zc].S() *lattice[xc][PBXR(PBYD(yc + Np, -lx, Np), 1, lx)][zc].S()) +
			(J1 / Ml) *(lattice[xc][PBYU(yc, lx, Np)][zc].S() + lattice[xc][PBYD(yc, -lx, Np)][zc].S()) +
			(J2 / Ml) *p +
			J3_prefac *(J3 / Ml) *(lattice[xc][yc][PBC(zc, 1, Nl)].S() + lattice[xc][yc][PBC(zc, -1, Nl)].S()) -
			K_prefac * Kx *(lattice[PBC(xc, 1, M)][yc][zc].S() + lattice[PBC(xc, -1, M)][yc][zc].S())
	);
	}
	else
	{
		dE = 2 *lattice[xc][yc][zc].S() *(-(J0 / Ml) *(lattice[xc][PBXR(yc, 1, lx)][zc].S() *lattice[xc][yc % Np][zc].S() *lattice[xc][PBYU(yc % Np, lx, Np)][zc].S() +
				lattice[xc][PBXL(yc, -1, lx)][zc].S() *lattice[xc][PBXL(yc % Np, -1, lx)][zc].S() *lattice[xc][PBXL(PBYU(yc % Np, lx, Np), -1, lx)][zc].S()) +
			(J1 / Ml) *(lattice[xc][PBXR(yc, 1, lx)][zc].S() + lattice[xc][PBXL(yc, -1, lx)][zc].S()) +
			(J2 / Ml) *p +
			J3_prefac *(J3 / Ml) *(lattice[xc][yc][PBC(zc, 1, Nl)].S() + lattice[xc][yc][PBC(zc, -1, Nl)].S()) -
			K_prefac * Kx *(lattice[PBC(xc, 1, M)][yc][zc].S() + lattice[PBC(xc, -1, M)][yc][zc].S())
	);
	}
	

	return dE;
}

deque<int> ITupdates::dip_enrg_diff_cluster(ITlattice ***lattice, ITlattice *plaquette, deque<int> cluster, deque<int> plaq_cluster, double *dE, int xc, int zc)
{
	double E_old, E_new, E_J3, px1, py1, px[4], py[4];
	int n[4];

	E_old = 0;
	for (std::deque<int>::iterator it = plaq_cluster.begin(); it != plaq_cluster.end(); ++it)
	{
		int yc = *it;
		plaquette[yc].set_S(0);
		n[0] = PBYU(yc, lx, Np);
		n[1] = PBXL(yc, -1, lx);
		n[2] = PBYD(yc, -lx, Np);
		n[3] = PBXR(yc, 1, lx);
		dip(lattice, xc, yc, zc, &px1, &py1);
		E_old += (J1 / Ml) *(lattice[xc][yc][zc].S() *lattice[xc][PBYU(yc, lx, Np)][zc].S() + 
						     lattice[xc][yc + Np][zc].S() *lattice[xc][PBXR(yc + Np, 1, lx)][zc].S());
		for (int i = 0; i < 4; ++i)
		{
			if (plaquette[n[i]].S() != 0)
			{
				dip(lattice, xc, n[i], zc, &px[i], &py[i]);
				E_old += (J2 / Ml) *(px[i] *px1 + py[i] *py1);
				E_old += (J1 / Ml) *(lattice[xc][n[i]][zc].S() *lattice[xc][PBYU(n[i], lx, Np)][zc].S() + 
						     lattice[xc][n[i] + Np][zc].S() *lattice[xc][PBXR(n[i] + Np, 1, lx)][zc].S());	
				plaquette[yc].reset_bond();
			}
		}
	}

		//if (myrank ==0) {cout << "   " << endl;} 				
	E_new = 0;	
	for (std::deque<int>::iterator it = plaq_cluster.begin(); it != plaq_cluster.end(); ++it)
	{
		int yc = *it;
		plaquette[yc].set_S(1);
		n[0] = PBYU(yc, lx, Np);
		n[1] = PBXL(yc, -1, lx);
		n[2] = PBYD(yc, -lx, Np);
		n[3] = PBXR(yc, 1, lx);
		dip(lattice, xc, yc, zc, &px1, &py1);
		E_new += (J1 / Ml) *(lattice[xc][yc][zc].S() *lattice[xc][PBYU(yc, lx, Np)][zc].S() + 
						     lattice[xc][yc + Np][zc].S() *lattice[xc][PBXR(yc + Np, 1, lx)][zc].S());		
		for (int i = 0; i < 4; ++i)
		{								
			if (plaquette[n[i]].S() != 1)
			{
				dip(lattice, xc, n[i], zc, &px[i], &py[i]);
				E_new += (J2 / Ml) *(px[i] *px1 + py[i] *py1);
				E_new += (J1 / Ml) *(lattice[xc][n[i]][zc].S() *lattice[xc][PBYU(n[i], lx, Np)][zc].S() + 
						     lattice[xc][n[i] + Np][zc].S() *lattice[xc][PBXR(n[i] + Np, 1, lx)][zc].S());
				//if (myrank ==0) {cout << " Ene is    " << *it << "  " << n[i] << "  " << plaquette[n[i]].S()<< endl;} 								     	
			}
		}
	}
	for (std::deque<int>::iterator it = plaq_cluster.begin(); it != plaq_cluster.end(); ++it)
	{
	plaquette[*it].set_S(-1);
	}	
	*dE = E_old - E_new;
				//if (myrank ==0) {cout << "       " << *dE << endl;} 								     	
	plaq_cluster.clear();
	return plaq_cluster;
}

void
ITupdates::pl(ITlattice ***lattice, int xc, int yc, int zc, double *pl_t)
{

	int l, k, dp[8];
	double Px[8], Py[8];
	if (yc < Np)
	{
		l = yc;
		k = PBYD(l, -lx, Np);

		dp[0] = PBXR(l, 1, lx);	
		dp[1] = PBYU(l, lx, Np);
		dp[2] = PBXL(l, -1, lx);	
		dp[3] = PBYD(l, -lx, Np);	
		dp[4] = PBXR(k, 1, lx);
		dp[5] = PBYU(k, lx, Np);
		dp[6] = PBXL(k, -1, lx);
		dp[7] = PBYD(k, -lx, Np);

		for (int b = 0; b < 8; b++) dip(lattice, xc, dp[b], zc, &Px[b], &Py[b]);
		*pl_t = 1 *(Px[0] + Px[1] + Px[2] + Px[3]) + 1 *(Py[0] + Py[1] + Py[2] + Py[3]) + -1 *(Px[4] + Px[5] + Px[6] + Px[7]) - 1 *(Py[4] + Py[5] + Py[6] + Py[7]);
		*pl_t /= 4;
	}		
	else
	{
		l = yc % Np;	
		k = PBXL(l, -1, lx);

		dp[0] = PBYD(l, -lx, Np);	
		dp[1] = PBXR(l, 1, lx);	
		dp[2] = PBYU(l, lx, Np);	
		dp[3] = PBXL(l, -1, lx);	
		dp[4] = PBYD(k, -lx, Np);
		dp[5] = PBXR(k, 1, lx);
		dp[6] = PBYU(k, lx, Np);
		dp[7] = PBXL(k, -1, lx);
		for (int b = 0; b < 8; b++) dip(lattice, xc, dp[b], zc, &Px[b], &Py[b]);
		*pl_t = 1 *(Px[0] + Px[1] + Px[2] + Px[3]) - 1 *(Py[0] + Py[1] + Py[2] + Py[3]) + -1 *(Px[4] + Px[5] + Px[6] + Px[7]) + 1 *(Py[4] + Py[5] + Py[6] + Py[7]);
		*pl_t /= 4;	       
	}
}



void
ITupdates::set_ferro_trial_state(ITlattice *** lattice){

  	for (int yc = 0; yc < Np; ++yc) {
  	  for (int xc = 0; xc < M; ++xc) {
  	    for (int zc = 0; zc < Nl; ++zc) { 
        	if (pow(-1, (yc % lx) + (yc / lx)) == 1){  
             		lattice[xc][yc][zc].set_S(-1);
			lattice[xc][PBYU(yc, lx, Np)][zc].set_S(1);
			lattice[xc][yc + Np][zc].set_S(-1);
			lattice[xc][PBXR(yc + Np, 1, lx)][zc].set_S(1);                        
        	}
  	    }
  	  }
  	}
}


void
ITupdates::dip(ITlattice ***lattice, int xc, int yc, int zc, double *px, double *py)
{
	*px = 0.25 * //pow(-1, (yc % lx) + (yc / lx)) *
		(lattice[xc][yc][zc].S() + lattice[xc][yc + Np][zc].S() - lattice[xc][PBXR(yc + Np, 1, lx)][zc].S() - lattice[xc][PBYU(yc, lx, Np)][zc].S());
	*py = 0.25 * //pow(-1, (yc % lx) + (yc / lx)) *
		(lattice[xc][yc][zc].S() + lattice[xc][PBXR(yc + Np, 1, lx)][zc].S() - lattice[xc][yc + Np][zc].S() - lattice[xc][PBYU(yc, lx, Np)][zc].S());
}
