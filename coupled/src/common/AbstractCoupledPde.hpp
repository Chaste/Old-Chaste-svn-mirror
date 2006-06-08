#ifndef ABSTRACTCOUPLEDPDE_HPP_
#define ABSTRACTCOUPLEDPDE_HPP_

#include "AbstractLinearParabolicPde.hpp"
#include <vector>
#include <petscvec.h>
#include "ReplicatableVector.hpp"
//#include <iostream>



template <int SPACE_DIM>
class AbstractCoupledPde : public AbstractLinearParabolicPde<SPACE_DIM>
{
protected:
    // timestep used by the pde solver
    double mBigTimeStep;

    // number of nodes in the mesh 
    int mNumNodes;
       
    // Lowest value of index that this part of the global object stores
    int mOwnershipRangeLo;
        
    // One more than the local highest index
    int mOwnershipRangeHi;
    
public:   
    AbstractCoupledPde(int numNodes, double bigTimeStep)
    {
        assert(numNodes > 0);
        
        mNumNodes = numNodes;
        mBigTimeStep = bigTimeStep;
        
        // Resize vectors to the appropriate size for each process:
        // Create a PETSc vector and use the ownership range of the PETSc vector
        // to size our C++ vectors
        Vec tempVec;
        VecCreate(PETSC_COMM_WORLD, &tempVec);
        VecSetSizes(tempVec, PETSC_DECIDE, numNodes);
        VecSetFromOptions(tempVec);
        VecGetOwnershipRange(tempVec,&mOwnershipRangeLo,&mOwnershipRangeHi);
        VecDestroy(tempVec); // vector no longer needed
    }
    
    void GetOwnershipRange(int &rLo, int &rHi)
    {
        rLo=mOwnershipRangeLo;
        rHi=mOwnershipRangeHi;
    }
        
};        
#endif /*ABSTRACTCOUPLEDPDE_HPP_*/
