#ifndef _NONLINEARHEATEQUATION4PDE_HPP_
#define _NONLINEARHEATEQUATION4PDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"


template <int SPACE_DIM>
class NonlinearHeatEquation4Pde : public AbstractNonlinearEllipticPde<SPACE_DIM>
{
public:

    double ComputeLinearSourceTerm(ChastePoint<SPACE_DIM> x)
    {
        return (-(exp(-2*x[0])-4*x[0]*exp(-2*x[0])+2*pow(x[0],2)*exp(-2*x[0])));
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<SPACE_DIM> , double )
    {
        return 0.0;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(ChastePoint<SPACE_DIM> , double u)
    {
        return identity_matrix<double>(SPACE_DIM) * u;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(ChastePoint<SPACE_DIM> , double )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
    
    double ComputeNonlinearSourceTermPrime(ChastePoint<SPACE_DIM> , double )
    {
        return 0.0;//(-(-4*exp(-x[0])+4*u));
    }
};

#endif //_NONLINEARHEATEQUATION4PDE_HPP_
