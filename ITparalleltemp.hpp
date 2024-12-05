#ifndef _ITPARALLELTEMPERING_HPP_DEFINED_
#define _ITPARALLELTEMPERING_HPP_DEFINED_
#include <iostream>
#include <algorithm>

using namespace std;
#pragma once
class ITparalleltemp {
  	public:
  	ITparalleltemp(){}
  	double Temp[96], Ene[96];
  	double d_Temp[96], dum_Temp[96];  
  	int ID[96];
	size_t v=1;
  		
      	void paralleltemp(int ,  ITupdates& );
      	void parallel_tempering(int ,  ITupdates& );
};
#endif
