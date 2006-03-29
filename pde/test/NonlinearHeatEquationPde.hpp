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

    double ComputeLinearSourceTerm(Point<SPACE_DIM> )
    {
    	return 1.0;
    }
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double )
    {
    	return 0.0;
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


#endif //_NONLINEARHEATEQUATIONPDE_HPP_
