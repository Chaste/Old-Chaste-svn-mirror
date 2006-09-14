#ifndef _ABSTRACTLINEARELLIPTICPDE_HPP_
#define _ABSTRACTLINEARELLIPTICPDE_HPP_

#include "UblasCustomFunctions.hpp"
#include "Point.hpp"
#include "Node.hpp"
#include <petscvec.h>


/**
 * AbstractLinearEllipticPde class.
 *
 * A general PDE of the form:
 * 0 = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 *
 * Parabolic PDEs are be derived from this (AbstractLinearParabolicPde)
 */

//// OLD NOTE: remember this if AbstractPde is brought back
// IMPORTANT NOTE: the inheritance of AbstractPde has to be 'virtual', ie 
// "class AbstractCardiacPde : public virtual AbstractPde"
// because AbstractPde will be the top class in a 'dreaded diamond':
//      A      
//     / \     A = AbstractPde, B = AbstractCardiac, 
//    B   C    C = AbtractLinearElliptic (and AbstractLinearParabolicPde) 
//     \ /     D = MonodomainPde
//      D
// 
// B and C must use virtual inheritence of A in order for D to only contain 1 instance
// of the member variables in A


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
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> x)=0;
    
    
    // - The following is defined in AbstractLinearParabolicPde, which inherits this class
    // - Compute the coefficient c(x) of du/dt
    //virtual double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> x)=0;
    
    virtual double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double u)
    {
        return ComputeNonlinearSourceTerm(node.GetPoint(), u);
    }
    
    virtual double ComputeLinearSourceTermAtNode(const Node<SPACE_DIM>& node)
    {
        return ComputeLinearSourceTerm(node.GetPoint());
    }
    
    virtual ~AbstractLinearEllipticPde()
    {
    
    }
};

#endif //_ABSTRACTLINEARELLIPTICPDE_HPP_
