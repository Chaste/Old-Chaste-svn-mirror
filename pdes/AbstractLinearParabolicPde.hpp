#ifndef _ABSTRACTLINEARPARABOLICPDE_HPP_
#define _ABSTRACTLINEARPARABOLICPDE_HPP_


/**
 * AbstractLinearParabolicPde class.
 * 
 * A general PDE of the form:
 * c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 * 
 */


#include "MatrixDouble.hpp"
#include "Point.hpp"
#include "Node.hpp"

template <int SPACE_DIM>
class AbstractLinearParabolicPde
{
public:
    virtual double ComputeLinearSourceTerm(Point<SPACE_DIM> x)=0;
    
    virtual double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x,
                                              double u)=0;

    virtual MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)=0;
    
    virtual double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> x)=0;
    
    virtual double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double u)
    {
        return ComputeNonlinearSourceTerm(node.GetPoint(), u);
    }
};
#endif //_ABSTRACTLINEARPARABOLICPDE_HPP_
