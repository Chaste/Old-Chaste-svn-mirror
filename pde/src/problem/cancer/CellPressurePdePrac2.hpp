#ifndef _CELLPRESSUREPDEPRAC2_HPP_
#define _CELLPRESSUREPDEPRAC2_HPP_

#include "AbstractLinearEllipticPde.hpp"

template <int SPACE_DIM>
class CellPressurePdePrac2:public AbstractLinearEllipticPde<SPACE_DIM>
{
	public:
    double mXalpha;
    double mXnecrotic;  /**< The boundary of the Necrotic Core */
    double mPressureAtBoundary;  /**< The pressure at the boundary of the necrotic and non-necrotic region */
    double mRho;
    
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
        
        double temp = 0.0;

        if(x[0] < mXalpha)
        {
            temp = -mRho;
        }
        else
        {
            temp = 1.0;
        }  
   
        return temp;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
        double temp = 0.0;
        if(x[0] < mXnecrotic)
        {
            temp = u - mPressureAtBoundary/2.0;
        }
        return temp;
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)
    {
    	return MatrixDouble::Identity(SPACE_DIM);
    }
};


#endif //_CELLPRESSUREPDEPRAC2_HPP_
