#ifndef TESTPARALLELVECTOR_HPP_
#define TESTPARALLELVECTOR_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

#include "ParallelVector.hpp"

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
        
        // initial condition;
        Vec striped;
        
        //unsigned big_lo, big_hi;;
        
        VecCreateMPI(PETSC_COMM_WORLD, 2*(hi-lo) , 2*10, &striped);
        double* p_striped;
        VecGetArray(striped, &p_striped);
        
        // stuff the vector with values: global index * local index
        double* p_vec;
        VecGetArray(vec, &p_vec);
        
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            p_vec[local_index] = local_index*global_index;
            p_striped[2*local_index  ] = local_index;
            p_striped[2*local_index+1] = global_index*global_index;
        }
        VecRestoreArray(vec, &p_vec);
        VecAssemblyBegin(vec);
        VecAssemblyEnd(vec);
        
        // create vec portion
        
        ParallelVector::SetProblemSize(10);

        ParallelVector parallel_vector(vec);
        ParallelVector parallel_striped(striped);
        ParallelVector::Stripe linear(parallel_striped,0);
        ParallelVector::Stripe quadratic(parallel_striped,1);
        
        // check that the portion has the correct range
        TS_ASSERT_EQUALS(ParallelVector::Begin().Global,lo);
        TS_ASSERT_EQUALS(ParallelVector::End().Global,hi);
        
        for (ParallelVector::Iterator index = ParallelVector::Begin();
             index!= ParallelVector::End();
             ++index)
        {
            TS_ASSERT_EQUALS(parallel_vector[index], index.Local*index.Global);
            TS_ASSERT_EQUALS(linear[index], index.Local);
            TS_ASSERT_EQUALS(quadratic[index], index.Global * index.Global);
        }
        
        ParallelVector::SetProblemSize(vec);
        
        // now write something to the vector - local_index * global_index
        for (ParallelVector::Iterator index = ParallelVector::Begin();
             index!= ParallelVector::End();
             ++index)
        {
            parallel_vector[index] =  -(double)(index.Local*index.Global);
            linear[index] =  1;
            quadratic[index] =  index.Local+1;
        }
        
        TS_ASSERT( ! (ParallelVector::Begin() == ParallelVector::End()) );        
        
        // ask the portion to restore the main vector
        parallel_vector.Restore();
        parallel_striped.Restore();
        
        // check that the vector was restored
        VecGetArray(vec, &p_vec);
        VecGetArray(striped, &p_striped);
        
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            TS_ASSERT_EQUALS(p_vec[local_index], -(double)local_index*global_index);
            TS_ASSERT_EQUALS(p_striped[2*local_index], (double)1);
            TS_ASSERT_EQUALS(p_striped[2*local_index+1], local_index+1);
        }

    }
           
        
};

#endif /*TESTPARALLELVECTOR_HPP_*/
