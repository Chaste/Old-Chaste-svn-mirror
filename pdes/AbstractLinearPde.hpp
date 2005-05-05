#ifndef _ABSTRACTLINEARPDE_HPP_
#define _ABSTRACTLINEARPDE_HPP_

#include "MatrixDouble.hpp"
#include "Point.hpp"
#include "Node.hpp"

/**
 * AbstractLinearPde class.
 * 
 * A general PDE of the form:
 * c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 * 
 * Both parabolic and elliptic PDEs can be derived from this.
 */

template <int SPACE_DIM>
class AbstractLinearPde
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
    
    /**
     * Compute the coefficient c(x) of du/dt
     */
	virtual double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> x)=0;


	virtual double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double u)
    {
        return ComputeNonlinearSourceTerm(node.GetPoint(), u);
    }
    
    virtual double ComputeLinearSourceTermAtNode(const Node<SPACE_DIM>& node)
    {
        return ComputeLinearSourceTerm(node.GetPoint());
    }

};


#endif //_ABSTRACTLINEARPDE_HPP_
