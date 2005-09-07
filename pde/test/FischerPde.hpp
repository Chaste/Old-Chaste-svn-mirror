#ifndef _FISCHERPDE_HPP_
#define _FISCHERPDE_HPP_


#include "AbstractLinearParabolicPde.hpp"
#include "Point.hpp"
#include <iostream>
/**
 * A simple parabolic PDE used in tests.
 */

template <int SPACE_DIM>
class FischerPde : public AbstractLinearParabolicPde<SPACE_DIM>
{
	
public:
	//std::vector<double> inputCache; 
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return 0.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
    	return u*(1-u);
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)
    {
    	return MatrixDouble::Identity(SPACE_DIM);
    }
    
	double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> x)
    {
    	return 1;
    }
    
    void PrepareForAssembleSystem(Vec v)
    {
        int lo, hi, numNodes;	
     	VecGetOwnershipRange(v, &lo, &hi);
        VecGetSize(v,&numNodes);
        double *vArray;
        VecGetArray(v, &vArray);
        double all_local_solutions[numNodes];
        
//        std::cout << "FischerPde::PrepareForAssembleSystem lo=" << lo
//            << " hi=" << hi << "numNodes="<< numNodes << std::endl;
        
        for (int i=0; i<numNodes; i++)
        {
        	if (lo <= i && i < hi)
	    	{ 
				all_local_solutions[i]=vArray[i-lo]; 
	        } 
	        else 
	        {
	           	all_local_solutions[i] =0.0;
	        }
        	
        }
 
    	double all_solutions[numNodes];
 		MPI_Allreduce(all_local_solutions, all_solutions, numNodes, MPI_DOUBLE, 
 		             MPI_SUM, PETSC_COMM_WORLD); 
    	
    		
    		AbstractLinearPde<SPACE_DIM>::inputCache.resize(numNodes);    
    	for (int i=0; i<numNodes; i++)
    	{
   			AbstractLinearPde<SPACE_DIM>::inputCache[i]=all_solutions[i];
    	}
    
     }
    
};

#endif
