#ifndef _EQUATIONFORNONNECROTICTUMOUR_HPP_
#define _EQUATIONFORNONNECROTICTUMOUR_HPP_

#include "AbstractNonlinearEllipticPde.hpp"
/*
 * this function is written for testing Practical 1 
 * 
 * 
 */

template <int SPACE_DIM>
class EquationForNonnecroticTumour : public AbstractNonlinearEllipticPde<SPACE_DIM>
{
public:

    double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
    {
        return 0.0;
    }
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
        const double sigma_0 = 1.0;
        const double sigma_1 = 1.0;
        
        return sigma_0*u/(1+sigma_1*u);
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x, double u)
    {        
        return MatrixDouble::Identity(SPACE_DIM) * 1.0;
    }
    
    MatrixDouble ComputeDiffusionTermPrime(Point<SPACE_DIM> x, double u)
    {
		return MatrixDouble::Identity(SPACE_DIM) * 0.0;
    }
    
    double ComputeNonlinearSourceTermPrime(Point<SPACE_DIM> x, double u)
    {
    	const double sigma_0 = 1.0;
        const double sigma_1 = 1.0;
        
    	return ((sigma_0*(1+sigma_1*u))-(sigma_1*sigma_0*u))/pow(1+sigma_1*u,2);
    }

};
#endif //_EQUATIONFORNONNECROTICTUMOUR_HPP_
