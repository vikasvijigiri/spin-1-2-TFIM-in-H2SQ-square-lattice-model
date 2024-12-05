#ifndef _ITUPDATES_HPP_DEFINED_
#define _ITUPDATES_HPP_DEFINED_
#include <deque>
#pragma once
using namespace std;
class ITupdates {
      	public: 
      	ITupdates(){}
      	
      	double E;
      	
  	int PBC(int, int, int);
  	int PBYU(int, int, int);
  	int PBYD(int, int, int);
  	int PBXR(int, int, int);
  	int PBXL(int, int, int);
  	
  	void pl(ITlattice *** , int, int, int, double * );
  	void dip(ITlattice *** , int, int, int, double * , double * );
  
	//Usual Imaginary-time QMC part
    	std::deque < int > dip_enrg_diff_cluster(ITlattice *** , ITlattice * , std::deque < int > , std::deque < int > , double * , int, int);
  	void make_bonds(ITlattice *** , ITlattice * , std::deque < int > * , std::deque < int > * , int, int, int, double);
  	double enrg_diff_local(ITlattice *** , int, int, int, int);
	double enrg(ITlattice *** ,  int );  	
  	int make_cluster(ITlattice *** , double);
  	void flip_cluster(ITlattice *** , std::deque < int > , bool, int, int, double );
  	void swendsen_flip(ITlattice *** , ITlattice * , double);


  	//Interlayer's part
  	int make_cluster_intl(ITlattice *** , int, double);
  	void make_bonds_intl(ITlattice *** , int, int, int, double);
  	void merge_clusters_intl(ITlattice *** , int, int);
  	void flip_cluster_intl(ITlattice *** , int, bool * );
  	double * cluster_field_intl(ITlattice *** , int, int);

};


inline int ITupdates::PBC(int a, int b, int L) { return (a + b + L) % L;}
inline int ITupdates::PBYU(int a, int b, int L) {return (a + b + L) % L;}
inline int ITupdates::PBYD(int a, int b, int L) {return (a % Np < lx) ? a + b + L : a + b;}
inline int ITupdates::PBXR(int a, int b, int L) {return ((a + b) % L == 0) ? a + b - L : a + b;}
inline int ITupdates::PBXL(int a, int b, int L) {return (a % L == 0) ? a + b + L : a + b;}
#endif



