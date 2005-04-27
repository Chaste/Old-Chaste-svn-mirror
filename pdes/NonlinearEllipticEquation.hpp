#ifndef _NONLINEARELLIPTICEQUATION_HPP_
#define _NONLINEARELLIPTICEQUATION_HPP_


#include "AbstractNonlinearEllipticPde.hpp"
/*
 * this function is written for testing Practical 1 
 * 
 * 
 */

template <int SPACE_DIM>
class NonlinearEllipticEquation : public AbstractNonlinearEllipticPde<SPACE_DIM>
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


#endif //_NONLINEARELLIPTICEQUATION_HPP_
