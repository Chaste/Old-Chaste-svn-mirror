#ifndef _TIMEDEPENDENTDIFFUSIONEQUATIONPDE_HPP_
#define _TIMEDEPENDENTDIFFUSIONEQUATIONPDE_HPP_

#include "AbstractLinearParabolicPde.hpp"
#include "ChastePoint.hpp"

/**
 * A simple parabolic PDE used in tests.
 */

template <int SPACE_DIM>
class TimeDependentDiffusionEquationPde : public AbstractLinearParabolicPde<SPACE_DIM>
{

public:
    double ComputeLinearSourceTerm(ChastePoint<SPACE_DIM> )
    {
        return 0.0;
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<SPACE_DIM> , double )
    {
        return 0.0;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(ChastePoint<SPACE_DIM> )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
    
    double ComputeDuDtCoefficientFunction(ChastePoint<SPACE_DIM> )
    {
        return 1;
    }
    
};

#endif //_TIMEDEPENDENTDIFFUSIONEQUATIONPDE_HPP_
