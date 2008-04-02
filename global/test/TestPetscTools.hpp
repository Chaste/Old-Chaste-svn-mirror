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

#ifndef TESTPETSCTOOLS_HPP_
#define TESTPETSCTOOLS_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestPetscTools : public CxxTest::TestSuite
{
public:
    void TestMostOfPetscTools()
    {
        PetscInt my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        TS_ASSERT_EQUALS(PetscTools::GetMyRank(), my_rank);
        bool am_master = (my_rank == 0);
        TS_ASSERT_EQUALS( PetscTools::AmMaster(), am_master);
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        bool is_sequential = (num_procs==1);
        TS_ASSERT_EQUALS( PetscTools::IsSequential(), is_sequential);

        ////////////////////////////////////////////////////
        // test CreateVec which returns a vec of constants
        ////////////////////////////////////////////////////
        Vec vec1 = PetscTools::CreateVec(10, 3.41);
        ReplicatableVector vec1_repl(vec1);
        
        TS_ASSERT_EQUALS(vec1_repl.size(), 10u);
        for(unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(vec1_repl[i], 3.41, 1e-12);
        }

        ////////////////////////////////////////////////////
        // test CreateVec which uses a std::vector of data
        ////////////////////////////////////////////////////
        std::vector<double> data(10);
        for(unsigned i=0; i<10; i++)
        {
            data[i] = i+0.45;
        }
        
        Vec vec2 = PetscTools::CreateVec(data);
        ReplicatableVector vec2_repl(vec2);
        
        TS_ASSERT_EQUALS(vec2_repl.size(), 10u);
        for(unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(vec2_repl[i], i+0.45, 1e-12);
        }
        
        ///////////////////////////////////////////////////
        // test SetupMatrix
        ///////////////////////////////////////////////////
        Mat mat;
        PetscTools::SetupMat(mat, 10, 11);
        int m,n;
        MatGetSize(mat, &m, &n);
        TS_ASSERT_EQUALS(m, 10);
        TS_ASSERT_EQUALS(n, 11);
        
        MatType type;
        MatGetType(mat,&type);
        //TS_ASSERT_EQUALS(type, MATMPIAIJ); // this does seem to work, but doesn't pass: it says "found (mpiaij != mpiaij)"

        VecDestroy(vec1);
        VecDestroy(vec2);
        MatDestroy(mat);
        
    }
    
    void TestBarrier()
    {
        // Testing the barrier method is kind of tricky, since we really want
        // to also check if it works when PETSc isn't set up.  So see TestPetscTools2.hpp!
        PetscTools::Barrier();
    }
    
    void TestReplicateError()
    {
        DistributedVector::SetProblemSize(1);
        if (DistributedVector::IsGlobalIndexLocal(0))
        {
            TS_ASSERT_THROWS_NOTHING(PetscTools::ReplicateException(true));
        }
        else
        {
            TS_ASSERT_THROWS_ANYTHING(PetscTools::ReplicateException(false));
        }
    }
};
#endif /*TESTPETSCTOOLS_HPP_*/
