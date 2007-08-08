#ifndef _ABSTRACTLINEARPARABOLICPDE_HPP_
#define _ABSTRACTLINEARPARABOLICPDE_HPP_

#include "AbstractLinearEllipticPde.hpp"
#include "ChastePoint.hpp"
#include "Node.hpp"

/**
 * AbstractLinearParabolicPde class.
 *
 * A general PDE of the form:
 * c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 *
 */
template <unsigned SPACE_DIM>
class AbstractLinearParabolicPde : public AbstractLinearEllipticPde<SPACE_DIM>
{
public:
    /**
     * The function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)" 
     */
    virtual double ComputeDuDtCoefficientFunction(ChastePoint<SPACE_DIM> x)=0;
    
    // The methods below are pure methods inherited from AbstractLinearEllipticPde:
    //virtual double ComputeLinearSourceTerm(ChastePoint<SPACE_DIM> x)=0;
    //virtual double ComputeNonlinearSourceTerm(ChastePoint<SPACE_DIM> x,
    //                                          double u)=0;
    //virtual MatrixDouble ComputeDiffusionTerm(ChastePoint<SPACE_DIM> x)=0;
    
    
};

#endif //_ABSTRACTLINEARPARABOLICPDE_HPP_
