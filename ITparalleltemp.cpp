#include "ITvariables.hpp"
#include "ITlattice.hpp"
#include "ITupdates.hpp"
#include "ITparalleltemp.hpp"
#include "mpi.h"
#include <iostream>
#include <algorithm>

static ran rann;
using namespace std;


void
ITparalleltemp::paralleltemp(int exc, ITupdates &upd)
{
        double ti,tR,tL,ei,eR,eL;
        int ii,iR,i,idi,idr,idlL, inc;



        //const auto ASIZE = std::extent<decltype(A)>::value;
	//const auto BSIZE = std::extent<decltype(B)>::value;
	const auto ASIZE = std::extent<decltype(Temp)>::value;
	const auto BSIZE = std::extent<decltype(Ene)>::value;  


        if (exc % 50 == 0) {
            MPI_Gather( & T, 1, MPI::DOUBLE, Temp, 1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
            MPI_Gather( & upd.E, 1, MPI::DOUBLE, Ene, 1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
            //MPI_Gather( & myrank, 1, MPI::INTEGER, ID, 1, MPI::INTEGER, 0, MPI::COMM_WORLD);
                //**************************************************************************************
                if (myrank == 0) {
                
                	/*
                	for (int i=0; i<numprocs; ++i){
			cout << " bef " << i << "   " << Temp[i] <<  "   "  << Ene[i] << endl;
			}
			*/
			// A and B must be the same size
			//assert(ASIZE == BSIZE);
					
			// p = sorted permutation of B largest to smallest

  			std::vector<size_t> p(ASIZE);	                             	
    			std::vector<size_t> pinv(ASIZE); 
			std::iota(p.begin(), p.end(), 0);
			std::sort(p.begin(), p.end(), [&](size_t i1, size_t i2) { return Ene[i1] < Ene[i2]; });
	              
			// pinv = inverse of p
 				
			for (size_t i = 0; i < ASIZE; ++i){
			    pinv[p[i]] = i;
			    //cout << pinv[p[i]] <<"  "  << p[i] << "   " <<Ene[pinv[p[i]]]<<  "  "<< Ene[p[i]] << endl;
			    }
			    
			
			 //Sort A largest to smallest
			 //std::sort(std::begin(Temp), std::end(Temp), [&](int v1, int v2) { return v1 > v2; }); 
			

			for (size_t i = 0; i < ASIZE; ++i){
			    dum_Temp[p[i]]=(pinv[i]+v)*0.02; 
			    //cout << pinv[p[i]] <<"  "  << p[i] << "   " <<Ene[pinv[p[i]]]<<  "  "<< Ene[p[i]] << "   "  << (pinv[i]+v)*0.01<< endl;
			    			    }
			for (int j=0; j<numprocs; ++j){
			  Temp[j]=dum_Temp[p[j]];
			  //cout << Temp[j] << "  " << Ene[j]<< endl;			
			}


			p.clear();
			pinv.clear();
			
			//p.shrink_to_fit();
			//pinv.shrink_to_fit();
			/*
			for (int i=0; i<numprocs; ++i){
			cout << " aft " << i << "   " << Temp[i] <<  "   "  << Ene[i] << endl;
			}
			*/
                } 
                // ***************************************!
                
                MPI_Scatter (Temp, 1, MPI::DOUBLE, & T,  1, MPI::DOUBLE, 0, MPI::COMM_WORLD); 
                //MPI_Scatter (ID, 1, MPI::INTEGER, & myrank, 1, MPI::INTEGER, 0, MPI::COMM_WORLD);   

        }
                
}


void
ITparalleltemp::parallel_tempering(int exc, ITupdates &upd) {

        double b1, b2, e1, e2;
        int c1, c2;
        int inc;

        if (exc % 10 == 0) {
                exc_index += 1;
           	MPI_Gather( & T, 1, MPI::DOUBLE, Temp, 1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
           	MPI_Gather( & upd.E, 1, MPI::DOUBLE, Ene, 1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
            	//MPI_Gather( & myrank, 1, MPI::INTEGER, ID, 1, MPI::INTEGER, 0, MPI::COMM_WORLD);
                //***************************************************************************************

                if (myrank == 0) {
                        for(int i=0; i<numprocs; ++i){ 
                        //int i = 5 - exc_index % 5;
                        inc=int(rann()*2) - 1;
                        if ((ID[i] + 1) == numprocs) inc = -1;  
			if ((ID[i]) == 0) inc = 1;                                             
			c1 = ID[i] + inc;
			c2 = ID[i];
                       
                        b1 = Temp[c1];
                        b2 = Temp[c2];
                        e1 = Ene[c1];
                        e2 = Ene[c2];


                        if (exp((1/b1 - 1/b2) * (e1 - e2)) > rann()) {
                                Temp[c2] = b1;
                                Temp[c1] = b2;
				ID[c1] = c2;
				ID[c2] = c1;
                        }
                        }
                }
                // ***************************************!
                MPI_Scatter (Temp, 1, MPI::DOUBLE, & T,  1, MPI::DOUBLE, 0, MPI::COMM_WORLD); 
                //MPI_Scatter (ID, 1, MPI::INTEGER, & myrank, 1, MPI::INTEGER, 0, MPI::COMM_WORLD);
        }
        // MPI_Barrier (MPI_COMM_WORLD, ierr);
        // ***************************************!
}
