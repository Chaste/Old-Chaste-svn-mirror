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

    SimpleNutrientPde(double coefficient);
    
    double ComputeConstantInUSourceTerm(const ChastePoint<DIM>& x);
    
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& x, Element<DIM,DIM>*);
    
    c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& );
    
};

template<unsigned DIM>
SimpleNutrientPde<DIM>::SimpleNutrientPde(double coefficient)
    : mCoefficient(coefficient)
{
}

template<unsigned DIM>
double SimpleNutrientPde<DIM>::ComputeConstantInUSourceTerm(const ChastePoint<DIM>& x)
{
    return 0.0;
}

template<unsigned DIM>
double SimpleNutrientPde<DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& x, Element<DIM,DIM>*)
{
    return -mCoefficient;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> SimpleNutrientPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& )
{
    return identity_matrix<double>(DIM);
} 

#endif /*SIMPLENUTRIENTPDE_HPP_*/
