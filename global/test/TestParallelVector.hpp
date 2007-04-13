#ifndef TESTPARALLELVECTOR_HPP_
#define TESTPARALLELVECTOR_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"
#include "ParallelProblem.hpp"
#include "ParallelVector.hpp"
//#include "VectorStripe.hpp"
#include "GlobalParallelProblem.hpp"

class TestParallelVector : public CxxTest::TestSuite
{
public:
    
    void TestReadAndRestore()
    {
        // create a 10 element petsc vector
        Vec vec;
        VecCreate(PETSC_COMM_WORLD, &vec);
        VecSetSizes(vec, PETSC_DECIDE, 10);
        VecSetFromOptions(vec);       
        // calculate my range
        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(vec,&petsc_lo,&petsc_hi);
        unsigned lo=(unsigned)petsc_lo;
        unsigned hi=(unsigned)petsc_hi;   
        // stuff the vector with values: global index * local index
        double* p_vec;
        VecGetArray(vec, &p_vec);
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            p_vec[local_index] = local_index*global_index;
        }
        VecRestoreArray(vec, &p_vec);
        VecAssemblyBegin(vec);
        VecAssemblyEnd(vec);
        //create a 20 element petsc vector with two stripes
        Vec striped;
        VecCreateMPI(PETSC_COMM_WORLD, 2*(hi-lo) , 2*10, &striped);
        // stuff the vectors
        double *p_striped;
        VecGetArray(striped, &p_striped);
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            p_striped[2*local_index] = global_index;
            p_striped[2*local_index+1] = -(double)global_index;
        }
        VecRestoreArray(striped, &p_striped);
        VecAssemblyBegin(striped);
        VecAssemblyEnd(striped);        
        
        
        
        //set the problem size
        gProblem.SetProblemSize(10);
        // create parallel vector
        ParallelVector vector(vec);
        ParallelVector big_vector(striped);
//        VectorStripe positive(big_vector,0);
//        VectorStripe negative(big_vector,1);
        // check we can read the values
        for (ParallelIterator index = gProblem.Begin();
             index!= gProblem.End();
             ++index)
        {
            TS_ASSERT_EQUALS(vector[index], gProblem.Local(index)*gProblem.Global(index));
        }
        // now write something to the vector: - local_index * global_index
        for (ParallelIterator index = gProblem.Begin();
             index!= gProblem.End();
             ++index)
        {
            vector[index] =  -(double)(gProblem.Local(index)*gProblem.Global(index));
        }       
        // restore the PETSc vector
        vector.Restore();
        
        // check that the PETSc vector was restored
        VecGetArray(vec, &p_vec);
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            TS_ASSERT_EQUALS(p_vec[local_index], -(double)local_index*global_index);
        }
    }
};

#endif /*TESTNEWVECTORPORTION_HPP_*/
