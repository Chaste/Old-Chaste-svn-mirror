#ifndef _ABSTRACTLINEARPARABOLICPDE_HPP_
#define _ABSTRACTLINEARPARABOLICPDE_HPP_

#include "AbstractLinearPde.hpp"
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
class AbstractLinearParabolicPde : public AbstractLinearPde<SPACE_DIM>
{
public:

	// The methods below are as for the base class, so commented out.
	//virtual double ComputeLinearSourceTerm(Point<SPACE_DIM> x)=0;
	//virtual double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x,
	//                                          double u)=0;
	//virtual MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)=0;
	//virtual double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> x)=0;
	//virtual double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double u)
	//virtual double ComputeLinearSourceTermAtNode(const Node<SPACE_DIM>& node)

    
};

#endif //_ABSTRACTLINEARPARABOLICPDE_HPP_
