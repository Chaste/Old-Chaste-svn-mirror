#ifndef _PRACTICALQUESTION2PDE_HPP_
#define _PRACTICALQUESTION2PDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"

template <int SPACE_DIM>
class Practical1Question2Pde:public AbstractNonlinearEllipticPde<SPACE_DIM>
{
	public:
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return 0.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
        double sigma_zero = 1.0;
        double sigma_one = 1.0;
    	return -(sigma_zero * u)/(1.0 + sigma_one * u);
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x, double u)
    {
    	return MatrixDouble::Identity(SPACE_DIM);
    }
};


#endif //_PRACTICALQUESTION2PDE_HPP_
