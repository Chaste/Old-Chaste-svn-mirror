#ifndef TESTREPLICATABLEVECTOR_HPP_
#define TESTREPLICATABLEVECTOR_HPP_
#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"

static const int VEC_SIZE=10;

class TestReplicatableVector : public CxxTest::TestSuite
{
public:
    void TestBasics()
    {
        ReplicatableVector rep_vector(VEC_SIZE);
        rep_vector[0]=15;
        rep_vector[1]=20;
        
        TS_ASSERT_EQUALS(rep_vector[0], 15);
        TS_ASSERT_EQUALS(rep_vector[1], 20);
        TS_ASSERT_EQUALS(rep_vector.size(), (unsigned) VEC_SIZE);
        
        rep_vector.resize(5);
        TS_ASSERT_EQUALS(rep_vector.size(), 5);
    }
    
    void TestReplication()
    {
        for (int vec_size=0; vec_size<10; vec_size++)
        {
            int lo, hi;
            Vec temp_vec;
            VecCreate(PETSC_COMM_WORLD, &temp_vec);
            VecSetSizes(temp_vec, PETSC_DECIDE, vec_size);
            VecSetFromOptions(temp_vec);
            VecGetOwnershipRange(temp_vec,&lo,&hi);
            VecDestroy(temp_vec); // vector no longer needed
            
            ReplicatableVector rep_vector(vec_size);
            for (int i=0; i<vec_size; i++)
            {
                rep_vector[i]=lo;
            }
            
            rep_vector.Replicate(lo, hi);
            
            for (int i=0; i<vec_size; i++)
            {
                if (lo<=i && i<hi)
                {
                    TS_ASSERT_EQUALS(rep_vector[i], lo);
                } else {
                    TS_ASSERT_DIFFERS(rep_vector[i], lo);
                }
            }
        }
    }
    
    void TestPetscReplication()
    {
        int lo, hi;
        Vec petsc_vec;
        VecCreate(PETSC_COMM_WORLD, &petsc_vec);
        VecSetSizes(petsc_vec, PETSC_DECIDE, VEC_SIZE);
        VecSetFromOptions(petsc_vec);
        VecGetOwnershipRange(petsc_vec,&lo,&hi);
        
        double *p_petsc_vec;
        
        VecGetArray(petsc_vec, &p_petsc_vec);
        for (int i=lo; i<hi; i++)
        {
            p_petsc_vec[i-lo]=lo;
        }
        VecRestoreArray(petsc_vec, &p_petsc_vec);
        VecAssemblyBegin(petsc_vec);
        VecAssemblyEnd(petsc_vec);
        
        ReplicatableVector rep_vec;
        rep_vec.ReplicatePetscVector(petsc_vec);    
            
        for (int i=0; i<VEC_SIZE; i++)
        {
            if (lo<=i && i<hi)
            {
                TS_ASSERT_EQUALS(rep_vec[i], lo);
            } else {
                TS_ASSERT_DIFFERS(rep_vec[i], lo);
            }
        }
        
        VecDestroy(petsc_vec);
    }
};
#endif /*TESTREPLICATABLEVECTOR_HPP_*/
