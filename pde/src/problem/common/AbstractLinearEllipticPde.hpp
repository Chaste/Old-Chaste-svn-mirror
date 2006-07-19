#ifndef _ABSTRACTLINEARELLIPTICPDE_HPP_
#define _ABSTRACTLINEARELLIPTICPDE_HPP_

#include "Point.hpp"
#include "AbstractLinearPde.hpp"

/**
 * AbstractLinearEllipticPde class.
 * 
 * A general PDE of the form:
 * 0 = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 * 
 */

template <int SPACE_DIM>
class AbstractLinearEllipticPde : public AbstractLinearPde<SPACE_DIM>
{
public:

	// The methods below are as for the base class, so commented out.
	//virtual double ComputeLinearSourceTerm(Point<SPACE_DIM> x)=0;
	//virtual double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x,
	//                                          double u)=0;
	//virtual MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)=0;

	/**
	 * An elliptic PDE has no time dependence, so the coefficient of du/dt is zero.
     * This method is never used.
	 */
    #define COVERAGE_IGNORE
	double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> )
	{
		return 0.0;
	}
    #undef COVERAGE_IGNORE
};

#endif //_ABSTRACTLINEARELLIPTICPDE_HPP_
