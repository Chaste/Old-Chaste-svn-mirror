#ifndef _ABSTRACTLINEARELLIPTICPDE_HPP_
#define _ABSTRACTLINEARELLIPTICPDE_HPP_

/**
 * AbstractLinearEllipticPde class.
 * 
 * A general PDE of the form:
 * 0 = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 * 
 */


#include "MatrixDouble.hpp"
#include "Point.hpp"

template <int SPACE_DIM>
class AbstractLinearEllipticPde
{
public:
/**
 * Compute Linear Source Term.
 * @param x The point in space at which the Linear Source Term is computed.
 */
    virtual double ComputeLinearSourceTerm(Point<SPACE_DIM> x)=0;

/**
 * Compute Nonlinear Source Term.
 * @param x The point in space at which the Nonlinear Source Term is computed.
 */
 
    virtual double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x,
                                              double u)=0;
/**
 * Compute Diffusion Term.
 * @param x The point in space at which the Diffusion Term is computed.
 * @return A matrix. 
 */
    virtual MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)=0;
};

#endif //_ABSTRACTLINEARELLIPTICPDE_HPP_
