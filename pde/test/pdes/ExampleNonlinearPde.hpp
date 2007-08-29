#ifndef _EXAMPLENONLINEARPDE_HPP_
#define _EXAMPLENONLINEARPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"
/**
 * An example 1D PDE: d/dx (u*du/dx) = -1
 */
class ExampleNonlinearPde:public AbstractNonlinearEllipticPde<1>
{

public:
    double ComputeLinearSourceTerm(ChastePoint<1> x)
    {
        return 1.0;
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<1> x, double u)
    {
        return 0.0;
    }
    
    c_matrix<double, 1, 1> ComputeDiffusionTerm(ChastePoint<1> , double u)
    {
        return identity_matrix<double>(1)*u;
    }
    
    c_matrix<double, 1, 1> ComputeDiffusionTermPrime(ChastePoint<1> , double )
    {
        return identity_matrix<double>(1)*0.0;
    }
    
    double ComputeNonlinearSourceTermPrime(ChastePoint<1> x, double u)
    {
        return 0.0;
    }
};

#endif //_EXAMPLENONLINEARPDE_HPP_
