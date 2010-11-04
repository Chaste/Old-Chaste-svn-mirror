/*
 * TestPetscVecTools.hpp
 *
 *  Created on: 4 Nov 2010
 *      Author: chaste
 */

#ifndef TESTPETSCVECTOOLS_HPP_
#define TESTPETSCVECTOOLS_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include "DistributedVectorFactory.hpp"
#include "ReplicatableVector.hpp"
#include "PetscVecTools.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedVector.hpp"

/**
 * Tests methods in the PETSc helper class PetscVecTools that are not directly related to LinearSystem
 */
class TestPetscVecTools : public CxxTest::TestSuite
{
public:

    void TestInterleavedVecScatter()
    {
        // Vectors will be twice PROBLEM_SIZE, since this is to be used in bidomain code.
        const unsigned PROBLEM_SIZE = 10;
        DistributedVectorFactory factory(PROBLEM_SIZE);

        // Source vector = [ 1 2 1 2 ... 1 2 ]
        Vec interleaved_vec=factory.CreateVec(2);
        DistributedVector interleaved_dist_vec = factory.CreateDistributedVector(interleaved_vec);
        DistributedVector::Stripe first_variable(interleaved_dist_vec, 0);
        DistributedVector::Stripe second_variable(interleaved_dist_vec, 1);
        for (DistributedVector::Iterator index = interleaved_dist_vec.Begin();
             index!= interleaved_dist_vec.End();
             ++index)
        {
            first_variable[index] = 1.0;
            second_variable[index] = 2.0;
        }

        // Destination vectors. It is important they have a compatible parallel layout with the source vector, hence the 2nd param of CreateVec()
        PetscInt local_num_rows;
        VecGetLocalSize(interleaved_vec, &local_num_rows);
        Vec first_variable_vec = PetscTools::CreateVec(PROBLEM_SIZE, local_num_rows/2);
        Vec second_variable_vec = PetscTools::CreateVec(PROBLEM_SIZE, local_num_rows/2);

        // Setup and perform scatter operation
        VecScatter first_variable_context;
        VecScatter second_variable_context;
        PetscVecTools::SetupInterleavedVectorScatterGather(interleaved_vec, first_variable_context, second_variable_context);
        PetscVecTools::DoInterleavedVecScatter(interleaved_vec, first_variable_context, first_variable_vec, second_variable_context, second_variable_vec);

        // Check destination vectors are made of 1s and 2s respectively.
        ReplicatableVector rep_1st_variable(first_variable_vec);
        ReplicatableVector rep_2nd_variable(second_variable_vec);
        TS_ASSERT_EQUALS(rep_1st_variable.GetSize(), PROBLEM_SIZE);
        TS_ASSERT_EQUALS(rep_2nd_variable.GetSize(), PROBLEM_SIZE);
        for (unsigned i=0; i<PROBLEM_SIZE; i++)
        {
            TS_ASSERT_EQUALS(rep_1st_variable[i], 1.0);
            TS_ASSERT_EQUALS(rep_2nd_variable[i], 2.0);
        }
    }

    void TestInterleavedVecGather()
    {
        // Vectors will be twice PROBLEM_SIZE, since this is to be used in bidomain code.
        const unsigned PROBLEM_SIZE = 10;
        DistributedVectorFactory factory(PROBLEM_SIZE);

        // Destination vector
        Vec interleaved_vec=factory.CreateVec(2);

        // Source vectors. It is important they have a compatible parallel layout with the destination vector, hence the 2nd param of CreateVec()
        PetscInt local_num_rows;
        VecGetLocalSize(interleaved_vec, &local_num_rows);
        Vec first_variable_vec = PetscTools::CreateVec(PROBLEM_SIZE, local_num_rows/2);
        Vec second_variable_vec = PetscTools::CreateVec(PROBLEM_SIZE, local_num_rows/2);

        // Fill in source vectors.
        DistributedVector dist_first_variable_vec = factory.CreateDistributedVector(first_variable_vec);
        DistributedVector dist_second_variable_vec = factory.CreateDistributedVector(second_variable_vec);
        for (DistributedVector::Iterator index = dist_first_variable_vec.Begin();
             index!= dist_first_variable_vec.End();
             ++index)
        {
            dist_first_variable_vec[index] = 1.0;
            dist_second_variable_vec[index] = 2.0;
        }

        // Setup and perform gather operation
        VecScatter first_variable_context;
        VecScatter second_variable_context;
        PetscVecTools::SetupInterleavedVectorScatterGather(interleaved_vec, first_variable_context, second_variable_context);
        PetscVecTools::DoInterleavedVecGather(interleaved_vec, first_variable_context, first_variable_vec, second_variable_context, second_variable_vec);

        // Check there are 1s and 2s interleaved in the destination vector
        ReplicatableVector rep_interleaved_vec(interleaved_vec);
        TS_ASSERT_EQUALS(rep_interleaved_vec.GetSize(), 2*PROBLEM_SIZE);
        for (unsigned i=0; i<PROBLEM_SIZE; i += 2)
        {
            TS_ASSERT_EQUALS(rep_interleaved_vec[i], 1);
            TS_ASSERT_EQUALS(rep_interleaved_vec[i+1], 2);
        }
    }
};

#endif /* TESTPETSCVECTOOLS_HPP_ */
