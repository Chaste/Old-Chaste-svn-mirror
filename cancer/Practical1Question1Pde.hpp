#ifndef _PRACTICALQUESTION1PDE_HPP_
#define _PRACTICALQUESTION1PDE_HPP_

#include "AbstractLinearEllipticPde.hpp"

template <int SPACE_DIM>
class Practical1Question1Pde:public AbstractLinearEllipticPde<SPACE_DIM>
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

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)
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


#endif //_PRACTICALQUESTION1PDE_HPP_
