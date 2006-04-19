#ifndef _CELLPRESSUREPDE_HPP_
#define _CELLPRESSUREPDE_HPP_

#include "AbstractLinearEllipticPde.hpp"

template <int SPACE_DIM>
class CellPressurePde:public AbstractLinearEllipticPde<SPACE_DIM>
{
private:    
    double mXalpha;
       
public:
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
        
        double temp = 0.0;
        double rho = 1.0;
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


#endif //_CELLPRESSUREPDE_HPP_
