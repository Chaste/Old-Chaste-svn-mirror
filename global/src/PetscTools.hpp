#ifndef PETSCTOOLS_HPP_
#define PETSCTOOLS_HPP_


#include <vector>
#include <petsc.h>
#include <petscvec.h>
#include <iostream>
#include <assert.h>
#include "DistributedVector.hpp"
#include "ReplicatableVector.hpp"

/**
 *  A helper class of static methods.
 */
class PetscTools
{
public :
    /**
     *  Just returns whether there is one process or not
     */
    static bool IsSequential()
    {
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        return (num_procs==1);
    }
 
    /**
     *  Create a vector of the specified size. SetFromOptions is called.
     */
    static Vec CreateVec(int size)
    {
        assert(size>0);
        Vec ret;
        VecCreate(PETSC_COMM_WORLD, &ret);
        VecSetSizes(ret, PETSC_DECIDE, size);
        VecSetFromOptions(ret);
        return ret;
    }
    
    /**
     *  Create a vector of the specified size with all values set to be the given
     *  constant. SetFromOptions is called.
     */
    static Vec CreateVec(int size, double value)
    {
        assert(size>0);
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
    
    /** 
     *  Create a Vec from the given data.
     */
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
    
    /**
     *  Set up a matrix - set the size using the given parameters, the type (default MATMPIAIJ). The
     *  number of local rows and columns is by default PETSC_DECIDE. SetFromOptions is called.
     */
    static void SetupMat(Mat& rMat, int numRows, int numColumns, MatType matType=MATMPIAIJ, int numLocalRows=PETSC_DECIDE, int numLocalColumns=PETSC_DECIDE)
    {
        assert(numRows>0);
        assert(numColumns>0);
        
        #if (PETSC_VERSION_MINOR == 2) //Old API
        MatCreate(PETSC_COMM_WORLD,numLocalRows,numLocalColumns,numRows,numColumns,&rMat);
        #else //New API
        MatCreate(PETSC_COMM_WORLD,&rMat);
        MatSetSizes(rMat,numLocalRows,numLocalColumns,numRows,numColumns);
        #endif
    
        MatSetType(rMat, matType);
        MatSetFromOptions(rMat);
    }
};
#endif /*PETSCTOOLS_HPP_*/
