#ifndef _TIMEDEPENDENTDIFFUSIONEQUATIONPDE_HPP_
#define _TIMEDEPENDENTDIFFUSIONEQUATIONPDE_HPP_

#include "AbstractLinearParabolicPde.hpp"
#include "Point.hpp"



template <int SPACE_DIM>
class TimeDependentDiffusionEquationPde : public AbstractLinearParabolicPde<SPACE_DIM>
{
	
public:
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return 0.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
    	return 0.0;
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

#endif //_TIMEDEPENDENTDIFFUSIONEQUATIONPDE_HPP_
