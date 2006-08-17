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
    unsigned mNumNodes;
    
    // Lowest value of index that this part of the global object stores
    unsigned mOwnershipRangeLo;
    
    // One more than the local highest index
    unsigned mOwnershipRangeHi;
    
public:
    AbstractCoupledPde(unsigned numNodes, double bigTimeStep)
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
        PetscInt temp_lo, temp_hi;
        VecGetOwnershipRange(tempVec, &temp_lo, &temp_hi);
        mOwnershipRangeLo=(unsigned) temp_lo;
        mOwnershipRangeHi=(unsigned) temp_hi;
        VecDestroy(tempVec); // vector no longer needed
    }
    
    void GetOwnershipRange(unsigned &rLo, unsigned &rHi)
    {
        rLo=mOwnershipRangeLo;
        rHi=mOwnershipRangeHi;
    }
    
};
#endif /*ABSTRACTCOUPLEDPDE_HPP_*/
