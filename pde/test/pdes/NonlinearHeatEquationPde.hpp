#ifndef _NONLINEARHEATEQUATIONPDE_HPP_
#define _NONLINEARHEATEQUATIONPDE_HPP_


#include "AbstractNonlinearEllipticPde.hpp"

/**
 *  A simple nonlinear elliptic PDE used by tests; Grad.(u Grad u) + 1 = 0
 */
template <int SPACE_DIM>
class NonlinearHeatEquationPde : public AbstractNonlinearEllipticPde<SPACE_DIM>
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
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& , double u)
    {
        return identity_matrix<double>(SPACE_DIM)*u;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(const ChastePoint<SPACE_DIM>& , double )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
    
    double ComputeNonlinearSourceTermPrime(const ChastePoint<SPACE_DIM>& , double )
    {
        return 0.0;
    }
    
    
};


#endif //_NONLINEARHEATEQUATIONPDE_HPP_
