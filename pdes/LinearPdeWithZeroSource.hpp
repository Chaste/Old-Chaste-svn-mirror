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
};

#endif //_LINEARPDEWITHZEROSOURCE_HPP_
