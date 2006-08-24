#ifndef _ABSTRACTLINEARPDE_HPP_
#define _ABSTRACTLINEARPDE_HPP_

#include "UblasCustomFunctions.hpp"
#include "AbstractPde.hpp"
#include "Point.hpp"
#include "Node.hpp"
#include <petscvec.h>



/**
 * AbstractLinearPde class.
 *
 * A general PDE of the form:
 * c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 *
 * Both parabolic and elliptic PDEs can be derived from this.
 */

// IMPORTANT NOTE: the inheritance of AbstractPde has to be 'virtual' because 
// AbstractPde will be the top class in a 'dreaded diamond':
//      A      
//     / \     A = AbstractPde, B = AbstractCardiac, C = AbtractLinearParabolic (etc)
//    B   C    D = MonodomainPde
//     \ /
//      D
// 
// B and C must use virtual inheritence of A in order for D to only contain 1 instance
// of the member variables in A
template <int SPACE_DIM>
class AbstractLinearPde : public virtual AbstractPde
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
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> x)=0;
    
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
    
    
    virtual ~AbstractLinearPde()
    {
    
    }
};


#endif //_ABSTRACTLINEARPDE_HPP_
