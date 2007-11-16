#ifndef _LINEARPDEWITHZEROSOURCE_HPP_
#define _LINEARPDEWITHZEROSOURCE_HPP_

#include "AbstractLinearEllipticPde.hpp"


/**
 * Linear PDE with zero source term and identity diffusion term.
 */
template <int SPACE_DIM>
class LinearPdeWithZeroSource:public AbstractLinearEllipticPde<SPACE_DIM>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& x)
    {
        return 0.0;
    }
    
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& )
    {
        return 0.0;
    }
        
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& x)
    {
        return identity_matrix<double>(SPACE_DIM);
    }
};

#endif //_LINEARPDEWITHZEROSOURCE_HPP_
