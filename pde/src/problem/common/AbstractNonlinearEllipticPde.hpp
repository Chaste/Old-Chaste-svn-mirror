#ifndef _ABSTRACTNONLINEARELLIPTICPDE_HPP_
#define _ABSTRACTNONLINEARELLIPTICPDE_HPP_

#include "UblasCustomFunctions.hpp"
#include "Point.hpp"


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

    virtual double ComputeLinearSourceTerm(Point<SPACE_DIM> x)=0;
    
    virtual double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x,
                                              double u)=0;
                                              
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> x,
            double u)=0;
            
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(Point<SPACE_DIM> x,
            double u)=0;
            
    virtual double ComputeNonlinearSourceTermPrime(Point<SPACE_DIM> x,
                                                   double u)=0;
    virtual ~AbstractNonlinearEllipticPde()
    {}
};

#endif //_ABSTRACTNONLINEARELLIPTICPDE_HPP_
