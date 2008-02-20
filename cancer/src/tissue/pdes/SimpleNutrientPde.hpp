#ifndef SIMPLENUTRIENTPDE_HPP_
#define SIMPLENUTRIENTPDE_HPP_

#include "AbstractLinearEllipticPde.hpp"

/**
 *  A simple nutrient PDE which is not directly coupled to the tissue.
 */
template<unsigned DIM>
class SimpleNutrientPde : public AbstractLinearEllipticPde<DIM>
{
private:
    double mCoefficient;

public:
    SimpleNutrientPde(double coefficient)
        : mCoefficient(coefficient)
    {
    }
    
    double ComputeConstantInUSourceTerm(const ChastePoint<DIM>& x)
    {
        return 0.0;
    }
    
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& x, Element<DIM,DIM>*)
    {
        return -mCoefficient;
    }
    
    c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& )
    {
        return identity_matrix<double>(DIM);
    }   
};

#endif /*SIMPLENUTRIENTPDE_HPP_*/
