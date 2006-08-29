#ifndef _ABSTRACTLINEARPARABOLICPDE_HPP_
#define _ABSTRACTLINEARPARABOLICPDE_HPP_

#include "AbstractLinearEllipticPde.hpp"
#include "Point.hpp"
#include "Node.hpp"

/**
 * AbstractLinearParabolicPde class.
 *
 * A general PDE of the form:
 * c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 *
 */
template <int SPACE_DIM>
class AbstractLinearParabolicPde : public AbstractLinearEllipticPde<SPACE_DIM>
{
public:
    /** 
     * The function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)" 
     */
    virtual double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> x)=0;

    // The methods below are pure methods inherited from AbstractLinearEllipticPde:
    //virtual double ComputeLinearSourceTerm(Point<SPACE_DIM> x)=0;
    //virtual double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x,
    //                                          double u)=0;
    //virtual MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)=0;

    
};

#endif //_ABSTRACTLINEARPARABOLICPDE_HPP_
