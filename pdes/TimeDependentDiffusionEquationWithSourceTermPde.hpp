#ifndef _TIMEDEPENDENTDIFFUSIONEQUATIONWITHSOURCETERMPDE_HPP_
#define _TIMEDEPENDENTDIFFUSIONEQUATIONWITHSOURCETERMPDE_HPP_


#include "AbstractLinearParabolicPde.hpp"
#include "Point.hpp"



template <int SPACE_DIM>
class TimeDependentDiffusionEquationWithSourceTermPde : public AbstractLinearParabolicPde<SPACE_DIM>
{
	
public:
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return 1.0;
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

#endif //_TIMEDEPENDENTDIFFUSIONEQUATIONWITHSOURCETERMPDE_HPP_
