#ifndef TESTDISTRIBUTEDVECTOR_HPP_
#define TESTDISTRIBUTEDVECTOR_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

#include "DistributedVector.hpp"

class TestDistributedVector : public CxxTest::TestSuite
{
public:
    
    void TestRead()
    {
        // WRITE VECTOR
        // create a 10 element petsc vector
        Vec vec;
        VecCreate(PETSC_COMM_WORLD, &vec);
        VecSetSizes(vec, PETSC_DECIDE, 10);
        VecSetFromOptions(vec);
        // calculate the range
        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(vec,&petsc_lo,&petsc_hi);
        unsigned lo=(unsigned)petsc_lo;
        unsigned hi=(unsigned)petsc_hi;   
        // create 20 element petsc vector
        Vec striped;
        VecCreateMPI(PETSC_COMM_WORLD, 2*(hi-lo) , 2*10, &striped);
        // write some values
        double* p_vec;
        VecGetArray(vec, &p_vec);
        double* p_striped;
        VecGetArray(striped, &p_striped);
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
        
        // READ VECTOR
        DistributedVector::SetProblemSize(10);
        DistributedVector distributed_vector(vec);
        DistributedVector distributed_vector2(striped);
        DistributedVector::Stripe linear(distributed_vector2,0);
        DistributedVector::Stripe quadratic(distributed_vector2,1);
        // check the range
        TS_ASSERT_EQUALS(DistributedVector::Begin().Global,lo);
        TS_ASSERT_EQUALS(DistributedVector::End().Global,hi);
        // read some values
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            TS_ASSERT_EQUALS(distributed_vector[index], index.Local*index.Global);
            TS_ASSERT_EQUALS(linear[index], index.Local);
            TS_ASSERT_EQUALS(quadratic[index], index.Global * index.Global);
        }

        
    }
    
    void TestWrite()
    {
        //WRITE VECTOR
        // create a 10 element petsc vector
        DistributedVector::SetProblemSize(10);
        Vec petsc_vec=DistributedVector::CreateVec();
        Vec striped=DistributedVector::CreateVec(2);
        DistributedVector distributed_vector(petsc_vec);
        DistributedVector distributed_vector2(striped);
        DistributedVector::Stripe linear(distributed_vector2,0);
        DistributedVector::Stripe quadratic(distributed_vector2,1);
        // write some values
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            distributed_vector[index] =  -(double)(index.Local*index.Global);
            linear[index] =  1;
            quadratic[index] =  index.Local+1;
        }
        
        distributed_vector.Restore();
        distributed_vector2.Restore();
        
        //READ VECTOR
        // calculate my range
        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(petsc_vec,&petsc_lo,&petsc_hi);
        unsigned lo=(unsigned)petsc_lo;
        unsigned hi=(unsigned)petsc_hi;   
        // read some values
        double* p_striped;
        VecGetArray(striped, &p_striped);
        double* p_vec;
        VecGetArray(petsc_vec, &p_vec);
        VecGetArray(petsc_vec, &p_vec);
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

#endif /*TESTDISTRIBUTEDVECTOR_HPP_*/
