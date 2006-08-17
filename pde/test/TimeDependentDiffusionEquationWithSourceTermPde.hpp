#ifndef _TIMEDEPENDENTDIFFUSIONEQUATIONWITHSOURCETERMPDE_HPP_
#define _TIMEDEPENDENTDIFFUSIONEQUATIONWITHSOURCETERMPDE_HPP_


#include "AbstractLinearParabolicPde.hpp"
#include "Point.hpp"

/**
 * A simple parabolic PDE used in tests.
 */

template <int SPACE_DIM>
class TimeDependentDiffusionEquationWithSourceTermPde : public AbstractLinearParabolicPde<SPACE_DIM>
{

public:
    double ComputeLinearSourceTerm(Point<SPACE_DIM> )
    {
        return 1.0;
    }
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double )
    {
        return 0.0;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
    
    double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> )
    {
        return 1;
    }
    
};

#endif //_TIMEDEPENDENTDIFFUSIONEQUATIONWITHSOURCETERMPDE_HPP_
