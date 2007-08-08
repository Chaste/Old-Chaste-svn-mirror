#ifndef _NONLINEARELLIPTICEQUATION_HPP_
#define _NONLINEARELLIPTICEQUATION_HPP_

#include "AbstractNonlinearEllipticPde.hpp"

/**
 *  A simple nonlinear elliptic PDE used by tests; Grad.(u Grad u) + 1 = 0
 */

template <int SPACE_DIM>
class NonlinearEllipticEquation : public AbstractNonlinearEllipticPde<SPACE_DIM>
{
public:

    double ComputeLinearSourceTerm(ChastePoint<SPACE_DIM> )
    {
        return 1.0;
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<SPACE_DIM> , double )
    {
        return 0.0;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(ChastePoint<SPACE_DIM> , double u)
    {
        return identity_matrix<double>(SPACE_DIM)*u;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(ChastePoint<SPACE_DIM> , double )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
    
    double ComputeNonlinearSourceTermPrime(ChastePoint<SPACE_DIM> , double )
    {
        return 0.0;
    }
    
};

#endif //_NONLINEARELLIPTICEQUATION_HPP_
