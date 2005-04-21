#ifndef _LINEARHEATEQUATIONPDE_HPP_
#define _LINEARHEATEQUATIONPDE_HPP_

#include "AbstractLinearEllipticPde.hpp"


template <int SPACE_DIM>
class LinearHeatEquationPde:public AbstractLinearEllipticPde<SPACE_DIM>
{
	public:
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return 1.0;
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


#endif //_HEATEQUATIONPDE_HPP_
