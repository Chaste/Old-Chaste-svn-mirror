#ifndef _ABSTRACTLINEARPDE_HPP_
#define _ABSTRACTLINEARPDE_HPP_

#include "MatrixDouble.hpp"
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

template <int SPACE_DIM>
class AbstractLinearPde
{
  public:
    // Kludge to make parallel stuff work...
    std::vector<double> inputCacheReplicated;
    

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
        ReplicateInputCache(currentSolution);
    }
 
    void ReplicateInputCache(Vec inputVector)
    {
        if (inputVector != NULL) 
        {
            int lo, hi, num_nodes;  
            VecGetOwnershipRange(inputVector, &lo, &hi);
            VecGetSize(inputVector,&num_nodes);       
            double *p_input_vector;
            VecGetArray(inputVector, &p_input_vector);
            double input_vector_local_array[num_nodes];
            for (int global_index=0; global_index<num_nodes; global_index++)
            {
                if (lo <= global_index && global_index < hi)
                { 
                	int local_index = global_index - lo;
                    input_vector_local_array[global_index]=p_input_vector[local_index]; 
                } 
                else 
                {
                    input_vector_local_array[global_index]=0.0;
                }
                
            }
     
            double input_vector_replicated_array[num_nodes];
            MPI_Allreduce(input_vector_local_array,input_vector_replicated_array, num_nodes, MPI_DOUBLE, 
                         MPI_SUM, PETSC_COMM_WORLD); 
            
            // Could be more efficient if MPI wrote to inputCacheReplicated above.
            inputCacheReplicated.resize(num_nodes);    
            for (int global_index=0; global_index<num_nodes; global_index++)
            {
                inputCacheReplicated[global_index]=input_vector_replicated_array[global_index];
            }
        } 
    }
    
    virtual ~AbstractLinearPde() 
    {
        
        
    }
};


#endif //_ABSTRACTLINEARPDE_HPP_
