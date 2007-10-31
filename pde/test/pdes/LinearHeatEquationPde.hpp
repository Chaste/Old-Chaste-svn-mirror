#ifndef _LINEARHEATEQUATIONPDE_HPP_
#define _LINEARHEATEQUATIONPDE_HPP_

#include "AbstractLinearEllipticPde.hpp"

/**
 * Steady state linear heat equation. Has unit source term and identity
 * diffusion term.
 */
template <int SPACE_DIM>
class LinearHeatEquationPde:public AbstractLinearEllipticPde<SPACE_DIM>
{
public:
    double ComputeLinearSourceTerm(ChastePoint<SPACE_DIM> )
    {
        return 1.0;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(ChastePoint<SPACE_DIM> )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
};

#endif //_HEATEQUATIONPDE_HPP_
