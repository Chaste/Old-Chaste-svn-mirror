#ifndef _ABSTRACTNONLINEARELLIPTICPDE_HPP_
#define _ABSTRACTNONLINEARELLIPTICPDE_HPP_

#include "UblasCustomFunctions.hpp"
#include "ChastePoint.hpp"


/**
 * AbstractNonlinearEllipticPde class.
 *
 * A simple elliptic PDE in 1 unknown with nonlinear diffusion term as
 * well as nonlinear source term:
 *
 *  0 = Grad.(DiffusionTerm(x,u)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 *
 */
template <unsigned SPACE_DIM>
class AbstractNonlinearEllipticPde
{
public:

    virtual double ComputeLinearSourceTerm(ChastePoint<SPACE_DIM> x)=0;
    
    virtual double ComputeNonlinearSourceTerm(ChastePoint<SPACE_DIM> x,
                                              double u)=0;
                                              
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(ChastePoint<SPACE_DIM> x,
            double u)=0;
            
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(ChastePoint<SPACE_DIM> x,
            double u)=0;
            
    virtual double ComputeNonlinearSourceTermPrime(ChastePoint<SPACE_DIM> x,
                                                   double u)=0;
    virtual ~AbstractNonlinearEllipticPde()
    {}
};

#endif //_ABSTRACTNONLINEARELLIPTICPDE_HPP_
