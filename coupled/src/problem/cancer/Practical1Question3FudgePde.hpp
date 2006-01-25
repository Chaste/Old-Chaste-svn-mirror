#ifndef _PRACTICALQUESTION3FUDGEPDE_HPP_
#define _PRACTICALQUESTION3FUDGEPDE_HPP_

#include "AbstractLinearEllipticPde.hpp"

template <int SPACE_DIM>
class Practical1Question3FudgePde:public AbstractLinearEllipticPde<SPACE_DIM>
{
	public:
    double mXalpha;
    double mXnew;
    
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
        double rho=1.0;
        double temp;
        if(x[0] < mXalpha)
        {
            temp = -(1.0/(mXnew*mXnew))*rho;
        }
        else
        {
            temp = (1.0/(mXnew*mXnew))*1.0;
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


#endif //_PRACTICALQUESTION3FUDGEPDE_HPP_
