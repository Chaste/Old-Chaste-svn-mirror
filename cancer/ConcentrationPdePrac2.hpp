#ifndef _CONCENTRATIONPDEPRAC2_HPP_
#define _CONCENTRATIONPDEPRAC2_HPP_

#include "AbstractNonlinearEllipticPde.hpp"
/*
 * Concentration Pde for Cancer Practicals One and Two
 *
*/
template <int SPACE_DIM>
class ConcentrationPdePrac2:public AbstractNonlinearEllipticPde<SPACE_DIM>
{
	private:
    // Function to calculate phi given an x value
    double CalculatePhi(Point<SPACE_DIM> x)
    {
        int lower = (int)((x[0]/mX)*mNumberOfNodes);
        double phi = mPhi[lower] + ( (x[0] - (mX/mNumberOfNodes)*lower) / 
                        (mX/mNumberOfNodes) ) * (mPhi[lower+1]-mPhi[lower]); 
        return phi;       
    }
    
    public:
    double mXnecrotic;  /**< The boundary of the Necrotic Core */
    double mX;          /**< The outer boundary of the Tumour */
    std::vector<double> mPhi;        /**< standard vector of the volume fraction of tumour cells. */
    int mNumberOfNodes;
    
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return 0.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
        double sigma_zero = 1.0;
        double sigma_one = 1.0;
        // Look at x, work out a value for phi from array pmPhi points to
        // Work out index of nodes around x
        double phi = CalculatePhi(x);
        double temp = -(sigma_zero * u)/(1.0 + sigma_one * u);
        if(x[0] < mXnecrotic)
        {
            temp = phi * temp ;
        }
    	return temp;
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
        double phi = CalculatePhi(x);
        double temp = (-(sigma_zero*(1+sigma_one*u))+(sigma_one*sigma_zero*u))/pow(1+sigma_one*u,2);
        if(x[0] < mXnecrotic)
        {
            temp = phi * temp ;
        }
    	return temp;
    }
};


#endif //_CONCENTRATIONPDEPRAC2_HPP_
