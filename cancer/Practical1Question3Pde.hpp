#ifndef _PRACTICALQUESTION3PDE_HPP_
#define _PRACTICALQUESTION3PDE_HPP_

#include "AbstractLinearEllipticPde.hpp"

template <int SPACE_DIM>
class Practical1Question3Pde:public AbstractLinearEllipticPde<SPACE_DIM>
{
	public:
    double mXalpha;
    
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
        double rho=1.0;
        double temp;
        if(x[0] < mXalpha)
        {
            temp = -rho;
        }
        else
        {
            temp = 1.0;
        }
		return temp;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
    	return 0.0;
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)
    {
    	return MatrixDouble::Identity(SPACE_DIM);
    }
};


#endif //_PRACTICALQUESTION3PDE_HPP_
