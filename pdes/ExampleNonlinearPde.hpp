#ifndef _EXAMPLENONLINEARPDE_HPP_
#define _EXAMPLENONLINEARPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"

/**
 * An example 1D PDE: d/dx (u*du/dx) = -1
 */
class ExampleNonlinearPde:public AbstractNonlinearEllipticPde<1>
{

public:
	double ComputeLinearSourceTerm(Point<1> x)
	{
		return 1.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<1> x, double u)
    {
    	return 0.0;
    }

    MatrixDouble ComputeDiffusionTerm(Point<1> x, double u)
    {
    	MatrixDouble I = MatrixDouble::Identity(1);
    	return u*I;
    }
    
    MatrixDouble ComputeDiffusionTermPrime(Point<1> x, double u)
    {
		return MatrixDouble::Identity(1) * 0.0;
    }
    
    double ComputeNonlinearSourceTermPrime(Point<1> x, double u)
    {
    	return 0.0;
    }
};

#endif //_EXAMPLENONLINEARPDE_HPP_
