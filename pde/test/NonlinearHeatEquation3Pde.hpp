#ifndef _NONLINEARHEATEQUATION3PDE_HPP_
#define _NONLINEARHEATEQUATION3PDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"

/**
 *  A simple nonlinear elliptic PDE used by tests; Grad.(u Grad u) - exp(-x) = 0
 */
 
template <int SPACE_DIM>
class NonlinearHeatEquation3Pde : public AbstractNonlinearEllipticPde<SPACE_DIM>
{
public:

    double ComputeLinearSourceTerm(Point<SPACE_DIM> )
    {
    	return 0.0;
    }
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double )
    {
    	return -exp(-x[0]);
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> , double u)
    {
		return MatrixDouble::Identity(SPACE_DIM) * u;
    }
    
    MatrixDouble ComputeDiffusionTermPrime(Point<SPACE_DIM> , double )
    {
		return MatrixDouble::Identity(SPACE_DIM) * 1.0;
    }
    
    double ComputeNonlinearSourceTermPrime(Point<SPACE_DIM> , double )
    {
    	return 0.0;
    }
};

#endif //_NONLINEARHEATEQUATION3PDE_HPP_
