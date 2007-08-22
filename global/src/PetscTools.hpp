#ifndef PETSCTOOLS_HPP_
#define PETSCTOOLS_HPP_


#include <vector>
#include <petscvec.h>
#include <iostream>
#include <assert.h>
#include "DistributedVector.hpp"
#include "ReplicatableVector.hpp"
#include "PetscSetupAndFinalize.hpp"

/**
 *  
 */
class PetscTools
{
public :
    static bool IsSequential()
    {
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        return (num_procs==1);
    }
 
    static Vec CreateVec(int size)
    {
        Vec ret;
        VecCreate(PETSC_COMM_WORLD, &ret);
        VecSetSizes(ret, PETSC_DECIDE, size);
        VecSetFromOptions(ret);
        return ret;
    }
    
    static Vec CreateVec(int size, double value)
    {
        Vec ret = CreateVec(size);
        
        #if (PETSC_VERSION_MINOR == 2) //Old API
        VecSet(&value, ret);
        #else
        VecSet(ret, value);
        #endif
        
        VecAssemblyBegin(ret);
        VecAssemblyEnd(ret);
        return ret;
    }
    
    static Vec CreateVec(std::vector<double> data)
    {
        assert(data.size()>0);
        Vec ret = CreateVec(data.size());

        double* p_ret;
        VecGetArray(ret, &p_ret);
        int lo, hi;
        VecGetOwnershipRange(ret, &lo, &hi);
        
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            p_ret[local_index] = data[global_index];
        }
        VecRestoreArray(ret, &p_ret);
        VecAssemblyBegin(ret);
        VecAssemblyEnd(ret);
        
        return ret;
    }
};
#endif /*PETSCTOOLS_HPP_*/
