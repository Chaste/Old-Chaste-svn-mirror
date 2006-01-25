#ifndef _CONCENTRATIONPDE_HPP_
#define _CONCENTRATIONPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"
/*
 * Concentration Pde for Cancer Practicals One and Two
 *
*/
template <int SPACE_DIM>
class ConcentrationPde:public AbstractNonlinearEllipticPde<SPACE_DIM>
{
	public:
    //double mXnecrotic;  /**< The boundary of the Necrotic Core */
    //double *pmPhi;        /**< The volume fraction of tumour cells. */
    
    
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
    
    MatrixDouble ComputeDiffusionTermPrime(Point<SPACE_DIM> x, double u)
    {
		return MatrixDouble::Identity(SPACE_DIM) * 0.0;
    }
    
    double ComputeNonlinearSourceTermPrime(Point<SPACE_DIM> x, double u)
    {
    	double sigma_zero = 1.0;
        double sigma_one = 1.0;
        return (-(sigma_zero*(1+sigma_one*u))+(sigma_one*sigma_zero*u))/pow(1+sigma_one*u,2);
    }
};


#endif //_CONCENTRATIONPDE_HPP_
