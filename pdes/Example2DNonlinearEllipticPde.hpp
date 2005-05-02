#ifndef _EXAMPLE2DNONLINEARELLIPTICPDE_HPP_
#define _EXAMPLE2DNONLINEARELLIPTICPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"
#include <cmath>

// template <int SPACE_DIM>
class Example2DNonlinearEllipticPde:public AbstractNonlinearEllipticPde<2>
{
	
	public:
	double ComputeLinearSourceTerm(Point<2> p)
	{
		double x = p[0], y = p[1];
    	return -( 4 + 8*x*x + 8*y*y );
	}
    
    double ComputeNonlinearSourceTerm(Point<2> p, double u)
    {
    	return 0.0;
    }

    MatrixDouble ComputeDiffusionTerm(Point<2> x, double u)
    {
    	MatrixDouble I = MatrixDouble::Identity(2);
    	return u*I;
    }
    
    MatrixDouble ComputeDiffusionTermPrime(Point<2> x, double u)
    {
		return MatrixDouble::Identity(2);
    }
    
    double ComputeNonlinearSourceTermPrime(Point<2> x, double u)
    {
    	return 0.0;
    }
};

#endif //_EXAMPLE2DNONLINEARELLIPTICPDE_HPP_
