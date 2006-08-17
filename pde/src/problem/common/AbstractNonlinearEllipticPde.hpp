#ifndef _ABSTRACTNONLINEARELLIPTICPDE_HPP_
#define _ABSTRACTNONLINEARELLIPTICPDE_HPP_

#include "UblasCustomFunctions.hpp"
#include "Point.hpp"

template <int SPACE_DIM>
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
