/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

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
        TS_ASSERT_EQUALS(rep_vector.size(), 5u);
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
            for (int global_index=0; global_index<vec_size; global_index++)
            {
                rep_vector[global_index]=lo;
            }
            
            rep_vector.Replicate(lo, hi);
            
            for (int global_index=0; global_index<vec_size; global_index++)
            {
                if (lo<=global_index && global_index<hi)
                {
                    TS_ASSERT_EQUALS(rep_vector[global_index], lo);
                }
                else
                {
                    TS_ASSERT_DIFFERS(rep_vector[global_index], lo);
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
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            p_petsc_vec[local_index]=lo;
        }
        VecRestoreArray(petsc_vec, &p_petsc_vec);
        VecAssemblyBegin(petsc_vec);
        VecAssemblyEnd(petsc_vec);
        
        ReplicatableVector rep_vec;
        rep_vec.ReplicatePetscVector(petsc_vec);
        
        for (int global_index=0; global_index<VEC_SIZE; global_index++)
        {
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_EQUALS(rep_vec[global_index], lo);
            }
            else
            {
                TS_ASSERT_DIFFERS(rep_vec[global_index], lo);
            }
        }
        
        VecDestroy(petsc_vec);
    }
    
    void TestPetscReplicationUsingAlternativeConstructor()
    {
        int lo, hi;
        Vec petsc_vec;
        VecCreate(PETSC_COMM_WORLD, &petsc_vec);
        VecSetSizes(petsc_vec, PETSC_DECIDE, VEC_SIZE);
        VecSetFromOptions(petsc_vec);
        VecGetOwnershipRange(petsc_vec,&lo,&hi);
        
        double *p_petsc_vec;
        
        VecGetArray(petsc_vec, &p_petsc_vec);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            p_petsc_vec[local_index]=lo;
        }
        VecRestoreArray(petsc_vec, &p_petsc_vec);
        VecAssemblyBegin(petsc_vec);
        VecAssemblyEnd(petsc_vec);
        
        ReplicatableVector rep_vec(petsc_vec);
        
        for (int global_index=0; global_index<VEC_SIZE; global_index++)
        {
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_EQUALS(rep_vec[global_index], lo);
            }
            else
            {
                TS_ASSERT_DIFFERS(rep_vec[global_index], lo);
            }
        }
        
        VecDestroy(petsc_vec);
    }
};
#endif /*TESTREPLICATABLEVECTOR_HPP_*/
