#ifndef _EXAMPLE2DNONLINEARELLIPTICPDE_HPP_
#define _EXAMPLE2DNONLINEARELLIPTICPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"

/**
 * A relatively simple PDE for testing the nonlinear elliptic assembler in 2D.
 */

class Example2DNonlinearEllipticPde:public AbstractNonlinearEllipticPde<2>
{

public:
    double ComputeLinearSourceTerm(ChastePoint<2> p)
    {
        double x = p[0], y = p[1];
        return -( 4 + 8*x*x + 8*y*y );
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<2> , double )
    {
        return 0.0;
    }
    
    c_matrix<double, 2, 2> ComputeDiffusionTerm(ChastePoint<2> , double u)
    {
        return identity_matrix<double>(2)*u;
    }
    
    c_matrix<double, 2, 2> ComputeDiffusionTermPrime(ChastePoint<2> , double )
    {
        return identity_matrix<double>(2);
    }
    
    double ComputeNonlinearSourceTermPrime(ChastePoint<2> , double )
    {
        return 0.0;
    }
    
    virtual ~Example2DNonlinearEllipticPde()
    {}
};

#endif //_EXAMPLE2DNONLINEARELLIPTICPDE_HPP_
