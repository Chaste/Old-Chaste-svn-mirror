#ifndef ELLIPTICPDEWITHLINEARSOURCE_HPP_
#define ELLIPTICPDEWITHLINEARSOURCE_HPP_

#include "AbstractLinearEllipticPde.hpp"

/**
 * The pde Div.(Grad u) + au + b = 0;
 */
template <int SPACE_DIM>
class EllipticPdeWithLinearSource :public AbstractLinearEllipticPde<SPACE_DIM>
{
private:
    double mCoeffOfU;
    double mConstant;

public:
    EllipticPdeWithLinearSource(double coeffOfU, double constant)
    {
        mCoeffOfU = coeffOfU;
        mConstant = constant;
    }

    double ComputeConstantInUSourceTerm(ChastePoint<SPACE_DIM> )
    {
        return mConstant;
    }

    double ComputeLinearInUCoeffInSourceTerm(ChastePoint<SPACE_DIM> )
    {
        return mCoeffOfU;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(ChastePoint<SPACE_DIM> )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
};


#endif /*ELLIPTICPDEWITHLINEARSOURCE_HPP_*/
