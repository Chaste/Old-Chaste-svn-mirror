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
	double ComputeLinearSourceTerm(Point<SPACE_DIM> )
	{
		return 0.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double u)
    {
    	return u*(1-u);
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
    	return MatrixDouble::Identity(SPACE_DIM);
    }
    
	double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> )
    {
    	return 1;
    }
    
  
};

#endif
