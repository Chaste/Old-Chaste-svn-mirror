#ifndef _NONLINEARHEATEQUATIONPDE_HPP_
#define _NONLINEARHEATEQUATIONPDE_HPP_


#include "AbstractNonlinearEllipticPde.hpp"


template <int SPACE_DIM>
class NonlinearHeatEquationPde : public AbstractNonlinearEllipticPde<SPACE_DIM>
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

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x, double u)
    {
		return MatrixDouble::Identity(SPACE_DIM) * u;
    }

};


#endif //_NONLINEARHEATEQUATIONPDE_HPP_
