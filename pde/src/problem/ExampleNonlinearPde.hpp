#ifndef _EXAMPLENONLINEARPDE_HPP_
#define _EXAMPLENONLINEARPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"
/**
 * An example 1D PDE: d/dx (u*du/dx) = -1
 *
 * This is here, rather than in pde/test/pdes, so that we have a 'real'
 * .cpp file in the pde component, in order that the pde library is
 * non-empty.  All the other .cpp files contain template classes, and
 * so are #included rather than linked against. 
 */
class ExampleNonlinearPde:public AbstractNonlinearEllipticPde<1>
{

public:
    double ComputeLinearSourceTerm(ChastePoint<1> x);
    
    double ComputeNonlinearSourceTerm(ChastePoint<1> x, double u);
    
    c_matrix<double, 1, 1> ComputeDiffusionTerm(ChastePoint<1> , double u);
    
    c_matrix<double, 1, 1> ComputeDiffusionTermPrime(ChastePoint<1> , double );
    
    double ComputeNonlinearSourceTermPrime(ChastePoint<1> x, double u);
};

#endif //_EXAMPLENONLINEARPDE_HPP_
