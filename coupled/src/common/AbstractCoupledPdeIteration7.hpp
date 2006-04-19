#ifndef ABSTRACTCOUPLEDPDEITERATION7_HPP_
#define ABSTRACTCOUPLEDPDEITERATION7_HPP_

#include "AbstractLinearParabolicPde.hpp"
#include <vector>
#include <petscvec.h>
//#include <iostream>

typedef std::vector<double> odeVariablesType;

template <int SPACE_DIM>
class AbstractCoupledPdeIteration7 : public AbstractLinearParabolicPde<SPACE_DIM>
{
protected:
    // timestep used by the pde solver
    double mBigTimeStep;

    // simulation time
    double mTime;   

    // number of nodes in the mesh 
    int mNumNodes;
       
    // Lowest value of index that this part of the global object stores
    int mOwnershipRangeLo;
        
    // One more than the local highest index
    int mOwnershipRangeHi;
    
    /** solutionCache stores the solutions to the ODEs (Icurrent) for
     *  each node in the global system.
     * 
     * This is replicated, i.e. use a global index for access.
     */
    std::vector<double> mSolutionCacheReplicated;
 
public:   
    //Constructor
    AbstractCoupledPdeIteration7(int numNodes, double tStart, double bigTimeStep)
    {
        assert(numNodes > 0);
        
        mNumNodes = numNodes;
        mBigTimeStep = bigTimeStep;
        mTime = tStart; /// \todo FIXME?: Is this overridden elsewhere?

        // Resize vectors to the appropriate size for each process:
        // Create a PETSc vector and use the ownership range of the PETSc vector
        // to size our C++ vectors
        Vec tempVec;
        VecCreate(PETSC_COMM_WORLD, &tempVec);
        VecSetSizes(tempVec, PETSC_DECIDE, numNodes);
        VecSetFromOptions(tempVec);
        VecGetOwnershipRange(tempVec,&mOwnershipRangeLo,&mOwnershipRangeHi);
        VecDestroy(tempVec); // vector no longer needed
        
        mSolutionCacheReplicated.resize(mNumNodes);
    }
 
    virtual void ReplicateSolutionCache(void)
    {
        double solution_cache_local_array[mNumNodes];
        for (int global_index=0; global_index<mNumNodes; global_index++)
        {
            if (mOwnershipRangeLo <= global_index && global_index < mOwnershipRangeHi)
            { 
                solution_cache_local_array[global_index]=mSolutionCacheReplicated[global_index];
            } 
            else 
            {
                solution_cache_local_array[global_index] =0.0;
            }
        }
 
        double solution_cache_replicated_array[mNumNodes];
        MPI_Allreduce(solution_cache_local_array, solution_cache_replicated_array, mNumNodes, MPI_DOUBLE, 
                     MPI_SUM, PETSC_COMM_WORLD); 
        
        // Could be more efficient if MPI wrote to solutionCacheReplicated above.
        for (int global_index=0; global_index<mNumNodes; global_index++)
        {
            mSolutionCacheReplicated[global_index]=solution_cache_replicated_array[global_index];
        }
    }
};        
#endif /*ABSTRACTCOUPLEDPDEITERATION7_HPP_*/
