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

    double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
    {
    	return 1.0;
    }
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
    	return 0.0;
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x, double u)
    {
		return MatrixDouble::Identity(SPACE_DIM);
    }

};


#endif //_NONLINEARLINEARHEATEQUATIONPDE_HPP_
