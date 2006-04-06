#ifndef ABSTRACTCOUPLEDPDEITERATION7_HPP_
#define ABSTRACTCOUPLEDPDEITERATION7_HPP_


#include "AbstractLinearParabolicPde.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include <vector>
#include <petscvec.h>
#include <iostream>

typedef std::vector<double> odeVariablesType;

template <int SPACE_DIM>
class AbstractCoupledPdeIteration7 : public AbstractLinearParabolicPde<SPACE_DIM>
{

public:
    // timestep used in the ode solvers        
 //   double mSmallTimeStep;

    // timestep used by the pde solver
    double mBigTimeStep;

    // simulation time
    double mTime;   

  //  AbstractIvpOdeSolver *mpOdeSolver;

    // number of nodes in the mesh 
    int mNumNodes;
        
    // Lowest value of index that this part of the global object stores
    int mOwnershipRangeLo;
        
    // One more than the local highest index
    int mOwnershipRangeHi;
    
    /** mOdeVarsAtNode[i] is a vector of the current values of the
     *  voltage, gating variables, intracellular Calcium concentration at node 
     *  i. The voltage returned by the ode solver is not used later since the pde
     *  solves for the voltage.
     * 
     * This is distributed, i.e. i should be a local index.
     */
//    std::vector<odeVariablesType>            mOdeVarsAtNode;
    
    
public:
    /** solutionCache stores the solutions to the ODEs (Icurrent) for
     *  each node in the global system.
     * 
     * This is replicated, i.e. use a global index for access.
     */
    std::vector<double> solutionCacheReplicated;
 
    // Replicated
    //std::vector<double>   inputCache;
 
    
    //Constructor
    AbstractCoupledPdeIteration7(int numNodes, /*AbstractIvpOdeSolver *pOdeSolver,*/ double tStart, double bigTimeStep/*, double smallTimeStep*/)
    {
 //       assert(smallTimeStep < bigTimeStep + 1e-10);
        assert(numNodes > 0);
        
        mNumNodes=numNodes;
        mBigTimeStep=bigTimeStep;
   //     mpOdeSolver=pOdeSolver;
   //     mSmallTimeStep=smallTimeStep;
     
        mTime = tStart;

        // Resize vectors to the appropriate size for each process:
        // Create a PETSc vector and use the ownership range of the PETSc vector
        // to size our C++ vectors
        Vec tempVec;
        VecCreate(PETSC_COMM_WORLD, &tempVec);
        VecSetSizes(tempVec, PETSC_DECIDE, numNodes);
        VecSetFromOptions(tempVec);
        VecGetOwnershipRange(tempVec,&mOwnershipRangeLo,&mOwnershipRangeHi);
        VecDestroy(tempVec); // no longer needed
        
//        mOdeVarsAtNode.resize(mOwnershipRangeHi-mOwnershipRangeLo);
      
        solutionCacheReplicated.resize(mNumNodes);


     }
     
 
     virtual void ReplicateSolutionCache(void)
     {
        
        double solution_cache_local_array[mNumNodes];
        for (int global_index=0; global_index<mNumNodes; global_index++)
        {
            if (mOwnershipRangeLo <= global_index && global_index < mOwnershipRangeHi)
            { 
                solution_cache_local_array[global_index]=solutionCacheReplicated[global_index];
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
            solutionCacheReplicated[global_index]=solution_cache_replicated_array[global_index];
        }
    
    }
    
//    odeVariablesType GetOdeVarsAtNode( int globalIndex )
//    {
//        if (!(mOwnershipRangeLo <= globalIndex && globalIndex < mOwnershipRangeHi)) {
//            std::cout << "i " << globalIndex << " lo " << mOwnershipRangeLo <<
 //               " hi " << mOwnershipRangeHi << std::endl;
 //       }
 //       assert(mOwnershipRangeLo <= globalIndex && globalIndex < mOwnershipRangeHi);
 //       return mOdeVarsAtNode[globalIndex-mOwnershipRangeLo];
 //   }
    
    /**
     * Apply same initial conditions to each node in the mesh
     */
//    void SetUniversalInitialConditions(odeVariablesType initialConditions)
  //  {
   //     for (int i=0; i<mOwnershipRangeHi-mOwnershipRangeLo; i++)
     //   {
       //     mOdeVarsAtNode[i] = initialConditions;
        //}
    //}

};        
#endif /*ABSTRACTCOUPLEDPDEITERATION7_HPP_*/
