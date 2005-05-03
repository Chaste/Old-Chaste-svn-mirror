#ifndef _PRACTICALQUESTION1PDENONLINEAR_HPP_
#define _PRACTICALQUESTION1PDENONLINEAR_HPP_

#include "AbstractNonlinearEllipticPde.hpp"

template <int SPACE_DIM>
class Practical1Question1PdeNonlinear:public AbstractNonlinearEllipticPde<SPACE_DIM>
{
	public:
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return -1.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
    	return 0.0;
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x, double u)
    {
    	return MatrixDouble::Identity(SPACE_DIM);
    }
    
    MatrixDouble ComputeDiffusionTermPrime(Point<SPACE_DIM> x, double u)
    {
		return MatrixDouble::Identity(SPACE_DIM) * 0.0;
    }
    
    double ComputeNonlinearSourceTermPrime(Point<SPACE_DIM> x, double u)
    {
    	return 0.0;
    }
};


#endif //_PRACTICALQUESTION1PDENONLINEAR_HPP_
