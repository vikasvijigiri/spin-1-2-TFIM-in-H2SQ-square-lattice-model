#ifndef _ITOBSERVABLES_HPP_DEFINED_
#define _ITOBSERVABLES_HPP_DEFINED_
#pragma once
class ITobservables {
  	public:
      	ITobservables(){}  	
        double op, op_sq, op_four, Ic, Ic_sq, E, E_sq;
      	double *data1, *data2, *wdata1, *wdata2;
      	int ** psite;
      	void observables(ITlattice *** ,  ITupdates & );
      	void binning_data(int , int );
      	void pbc();
      	void Initiate_observables();
};
#endif
