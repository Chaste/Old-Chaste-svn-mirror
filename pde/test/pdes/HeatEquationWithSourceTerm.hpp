#ifndef _HEATEQUATIONWITHSOURCETERM_HPP_
#define _HEATEQUATIONWITHSOURCETERM_HPP_


#include "AbstractLinearParabolicPde.hpp"
#include "ChastePoint.hpp"

/**
 * A simple parabolic PDE used in tests.
 */

template <int SPACE_DIM>
class HeatEquationWithSourceTerm : public AbstractLinearParabolicPde<SPACE_DIM>
{

public:
    double ComputeLinearSourceTerm(const ChastePoint<SPACE_DIM>& )
    {
        return 1.0;
    }
    
    double ComputeNonlinearSourceTerm(const ChastePoint<SPACE_DIM>& , double )
    {
        return 0.0;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& , Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
    {
        return identity_matrix<double>(SPACE_DIM);
    }
    
    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& )
    {
        return 1;
    }
    
};

#endif //_HEATEQUATIONWITHSOURCETERM_HPP_
