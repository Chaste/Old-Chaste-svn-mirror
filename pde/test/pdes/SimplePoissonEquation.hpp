#ifndef _SIMPLEPOISSONEQUATION_HPP_
#define _SIMPLEPOISSONEQUATION_HPP_

#include "AbstractLinearEllipticPde.hpp"

/**
 * Steady state linear heat equation. Has unit source term and identity
 * diffusion term.
 */
template <int SPACE_DIM>
class SimplePoissonEquation:public AbstractLinearEllipticPde<SPACE_DIM>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& )
    {
        return 1.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& )
    {
        return 0.0;
    }    

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
};

#endif //_SIMPLEPOISSONEQUATION_HPP_
