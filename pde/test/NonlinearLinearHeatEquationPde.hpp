#ifndef _NONLINEARLINEARHEATEQUATIONPDE_HPP_
#define _NONLINEARLINEARHEATEQUATIONPDE_HPP_


#include "AbstractNonlinearEllipticPde.hpp"

/**
 * Linear heat equation, using the nonlinear class format.
 */
template <int SPACE_DIM>
class NonlinearLinearHeatEquationPde : public AbstractNonlinearEllipticPde<SPACE_DIM>
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

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> , double )
    {
		return MatrixDouble::Identity(SPACE_DIM);
    }
    
    MatrixDouble ComputeDiffusionTermPrime(Point<SPACE_DIM> , double )
    {
		return MatrixDouble::Identity(SPACE_DIM) * 0.0;
    }
    
    double ComputeNonlinearSourceTermPrime(Point<SPACE_DIM> , double )
    {
    	return 0.0;
    }

};


#endif //_NONLINEARLINEARHEATEQUATIONPDE_HPP_
