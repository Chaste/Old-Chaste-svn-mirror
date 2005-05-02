#ifndef _NONLINEARHEATEQUATION2PDE_HPP_
#define _NONLINEARHEATEQUATION2PDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"


template <int SPACE_DIM>
class NonlinearHeatEquation2Pde : public AbstractNonlinearEllipticPde<SPACE_DIM>
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
		return MatrixDouble::Identity(SPACE_DIM) * (1.0/u);
    }
    
    MatrixDouble ComputeDiffusionTermPrime(Point<SPACE_DIM> x, double u)
    {
		return MatrixDouble::Identity(SPACE_DIM) * (-1.0/pow(u,2));
    }
    
    double ComputeNonlinearSourceTermPrime(Point<SPACE_DIM> x, double u)
    {
    	return 0.0;
    }
};

#endif //_NONLINEARHEATEQUATION2PDE_HPP_
