#ifndef _NONLINEARHEATEQUATION4PDE_HPP_
#define _NONLINEARHEATEQUATION4PDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"


template <int SPACE_DIM>
class NonlinearHeatEquation4Pde : public AbstractNonlinearEllipticPde<SPACE_DIM>
{
public:

    double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
    {
    	return (-(exp(-2*x[0])-4*x[0]*exp(-2*x[0])+2*pow(x[0],2)*exp(-2*x[0])));
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
    	return 0.0;//(-(-4*exp(-x[0])+4*u));
    }
};

#endif //_NONLINEARHEATEQUATION4PDE_HPP_
