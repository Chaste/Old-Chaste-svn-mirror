#ifndef _ABSTRACTLINEARPDE_HPP_
#define _ABSTRACTLINEARPDE_HPP_

#include "MatrixDouble.hpp"
#include "Point.hpp"
#include "Node.hpp"
#include "petscvec.h"

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
    // Kludge to make parallel stuff work...
    std::vector<double> inputCache;
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
    
    
    
    /**
    * PrepareForAssembleSystem is a virtual method.
    * It's called by the AssembleSystem method of the assembler before any other
    * useful work happens.  The idea is that a *coupled system* will want to 
    * solve all the ODE systems before the PDE is solved.  A *parallel* coupled
    * system will want to solve the ODE systems and distribute the answers 
    * before anything else happens.
    */ 
    virtual void PrepareForAssembleSystem(Vec currentSolution)
    {
        DistributeInputCache(currentSolution);
    }
 
    void DistributeInputCache(Vec inputVector)
    {
        if (inputVector != NULL) 
        {
            int lo, hi, num_nodes;  
            VecGetOwnershipRange(inputVector, &lo, &hi);
            VecGetSize(inputVector,&num_nodes);       
            double *vArray;
            VecGetArray(inputVector, &vArray);
            double all_local_solutions[num_nodes];
            for (int i=0; i<num_nodes; i++)
            {
                if (lo <= i && i < hi)
                { 
                    all_local_solutions[i]=vArray[i-lo]; 
                } 
                else 
                {
                    all_local_solutions[i]=0.0;
                }
                
            }
     
            double all_solutions[num_nodes];
            MPI_Allreduce(all_local_solutions, all_solutions, num_nodes, MPI_DOUBLE, 
                         MPI_SUM, PETSC_COMM_WORLD); 
            
            // Could be more efficient if MPI wrote to inputCache above.
            inputCache.resize(num_nodes);    
            for (int i=0; i<num_nodes; i++)
            {
                inputCache[i]=all_solutions[i];
            }
        } 
    }
};


#endif //_ABSTRACTLINEARPDE_HPP_
