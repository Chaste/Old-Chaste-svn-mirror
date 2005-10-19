#ifndef _FAKELR91PDE_HPP_
#define _FAKELR91PDE_HPP_


#include "AbstractLinearParabolicPde.hpp"
#include "Point.hpp"

/**
 * A simple parabolic PDE used in tests.
 */

template <int SPACE_DIM>
class FakeLr91Pde : public AbstractLinearParabolicPde<SPACE_DIM>
{
	
public:
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return 0.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
        if (x[0] == 0.0)
        {
            //return 598.591;
            return 0.0;
        }
        else
        {
    	    //return -0.0144045;
            return 0.0;
        }
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)
    {
    	return MatrixDouble::Identity(SPACE_DIM);
    }
    
	double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> x)
    {
    	return 1;
    }
    
};

#endif 
