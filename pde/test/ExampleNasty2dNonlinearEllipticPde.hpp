#ifndef _EXAMPLENASTY2DNONLINEARELLIPTICPDE_HPP_
#define _EXAMPLENASTY2DNONLINEARELLIPTICPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"
#include <cmath>

/**
 * A fairly nasty PDE for testing the nonlinear elliptic assembler in 2D.
 */

class ExampleNasty2dNonlinearEllipticPde:public AbstractNonlinearEllipticPde<2>
{

public:
    double ComputeLinearSourceTerm(ChastePoint<2> )
    {
        return 0;
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<2> p, double u)
    {
        double x = p[0], y = p[1];
        return -4*(u*cos(x)*cos(x) + sin(x)*sin(x)*cos(x)*cos(x) + y*y);
    }
    
    c_matrix<double, 2, 2> ComputeDiffusionTerm(ChastePoint<2> , double u)
    {
        return identity_matrix<double>(2)*u;
    }
    
    c_matrix<double, 2, 2> ComputeDiffusionTermPrime(ChastePoint<2> , double )
    {
        return identity_matrix<double>(2);
    }
    
    double ComputeNonlinearSourceTermPrime(ChastePoint<2> p, double )
    {
        return -(cos(p[0])*cos(p[0]));
    }
    
    virtual ~ExampleNasty2dNonlinearEllipticPde()
    {}
};

#endif //_EXAMPLENASTY2DNONLINEARELLIPTICPDE_HPP_
