#ifndef _PRACTICALQUESTION2FUDGEPDE_HPP_
#define _PRACTICALQUESTION2FUDGEPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"

template <int SPACE_DIM>
class Practical1Question2FudgePde:public AbstractNonlinearEllipticPde<SPACE_DIM>
{
	public:
    
    double mXnew;
    
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return 0.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
        double sigma_zero = 1.0;
        double sigma_one = 1.0;
        
    	return -(1.0/(mXnew*mXnew))*(sigma_zero * u)/(1.0 + sigma_one * u);
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
    	double sigma_zero = 1.0;
        double sigma_one = 1.0;
        
    	return (1.0/(mXnew*mXnew))*(-(sigma_zero*(1+sigma_one*u))+(sigma_one*sigma_zero*u))/pow(1+sigma_one*u,2);
    }
};


#endif //_PRACTICALQUESTION2FUDGEPDE_HPP_
