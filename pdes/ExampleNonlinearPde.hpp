#ifndef _EXAMPLENONLINEARPDE_HPP_
#define _EXAMPLENONLINEARPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"
//#include <cmath>

// template <int SPACE_DIM>
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
};

#endif //_EXAMPLENONLINEARPDE_HPP_
