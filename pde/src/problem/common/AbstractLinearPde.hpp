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
private:
    // Kludge to make parallel stuff work...
    std::vector<double> mInputCacheReplicated;
    
protected:

    /**
     * Replicate a vector over all processes.
     * 
     * Each process knows its local part of the vector.  This method shares that knowledge
     * across all the processes.
     * 
     * @param lo  The start of our ownership range
     * @param hi  One past the end of our ownership range
     * @param size  The size of the vector to be replicated
     * @param input_array  The local portion of the array to be replicated.  Should
     *    contain hi-lo entries.  (If your input vector is larger, just pass the address
     *    of the first entry in your ownership range.)
     * @param output_array  The array to store replicated data in.  May be the same as
     *    input_array.  Memory for size entries should have been allocated already.
     */
    void ReplicateVector(unsigned lo, unsigned hi, unsigned size,
                         double *input_array,
                         double *output_array)
    {
        // Set up an array for MPI replication to use
        double input_vector_local_array[size];
        for (unsigned global_index=0; global_index<size; global_index++)
        {
            if (lo <= global_index && global_index < hi)
            { 
                unsigned local_index = global_index - lo;
                input_vector_local_array[global_index] = input_array[local_index]; 
            } 
            else 
            {
                input_vector_local_array[global_index] = 0.0;
            }
        }
        
        // Replicate
        MPI_Allreduce(input_vector_local_array, output_array, size,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    }
    
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
        if (currentSolution != NULL) 
        {
            int lo, hi, num_nodes;  
            VecGetOwnershipRange(currentSolution, &lo, &hi);
            VecGetSize(currentSolution, &num_nodes);       
            double *input_array;
            VecGetArray(currentSolution, &input_array);
            
            mInputCacheReplicated.resize(num_nodes);
            ReplicateVector(lo, hi, num_nodes, input_array, &(mInputCacheReplicated[0]));
            
            VecRestoreArray(currentSolution, &input_array);
        } 
    }
    
    double GetInputCacheMember(unsigned int i)
    {
        assert(i<mInputCacheReplicated.size());
        return(mInputCacheReplicated[i]);
    }
    virtual ~AbstractLinearPde() 
    {
        
        
    }
};


#endif //_ABSTRACTLINEARPDE_HPP_
