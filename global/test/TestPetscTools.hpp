/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef TESTPETSCTOOLS_HPP_
#define TESTPETSCTOOLS_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>
#include "DistributedVectorFactory.hpp"
#include "ReplicatableVector.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"

class TestPetscTools : public CxxTest::TestSuite
{
public:
    void TestMostOfPetscTools()
    {
        PetscInt my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        TS_ASSERT_EQUALS(PetscTools::GetMyRank(), (unsigned)my_rank);
        bool am_master = (my_rank == 0);
        TS_ASSERT_EQUALS( PetscTools::AmMaster(), am_master);

        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        TS_ASSERT_EQUALS( PetscTools::GetNumProcs(), (unsigned)num_procs);
        bool is_sequential = (num_procs==1);
        TS_ASSERT_EQUALS( PetscTools::IsSequential(), is_sequential);

        ////////////////////////////////////////////////////
        // test CreateVec which returns a vec of constants
        ////////////////////////////////////////////////////
        Vec vec1 = PetscTools::CreateVec(10, 3.41);
        ReplicatableVector vec1_repl(vec1);

        TS_ASSERT_EQUALS(vec1_repl.size(), 10u);
        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(vec1_repl[i], 3.41, 1e-12);
        }

        ////////////////////////////////////////////////////
        // test CreateVec which uses a std::vector of data
        ////////////////////////////////////////////////////
        std::vector<double> data(10);
        for (unsigned i=0; i<10; i++)
        {
            data[i] = i+0.45;
        }

        Vec vec2 = PetscTools::CreateVec(data);

        ReplicatableVector vec2_repl(vec2);

        TS_ASSERT_EQUALS(vec2_repl.size(), 10u);
        for (unsigned i=0; i<10; i++)
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

#if (PETSC_VERSION_MAJOR == 3)
        const MatType type;
#else
        MatType type;
#endif
        MatGetType(mat,&type);
        TS_ASSERT(strcmp(type, MATMPIAIJ)==0);

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
        //
        // A factory is created here to set static members on DistributedVector.
        // As part of #988 these members are being moved to non-static members on
        // DistributedVectorFactory. This should be refactored once this is complete.
        //
        DistributedVectorFactory factory(1);
        if (factory.IsGlobalIndexLocal(0))
        {
            TS_ASSERT_THROWS_NOTHING(PetscTools::ReplicateException(true));
        }
        else
        {
            TS_ASSERT_THROWS_THIS(PetscTools::ReplicateException(false), "Another process threw an exception; bailing out.");
        }
    }

    void TestDumpPetscObjects()
    {
        Mat matrix;
        Vec vector;

        PetscTools::SetupMat(matrix, 10, 10, (MatType)MATMPIAIJ);

        VecCreate(PETSC_COMM_WORLD, &vector);
        VecSetSizes(vector, PETSC_DECIDE, 10);
        VecSetFromOptions(vector);

        PetscInt lo, hi;
        VecGetOwnershipRange(vector, &lo, &hi);

        for (int row=0; row<10; row++)
        {
            if (row >= lo && row < hi)
            {
                for (int col=0; col<10; col++)
                {
                    MatSetValue(matrix, row, col, (double) 10*row+col+1, INSERT_VALUES);
                }

                double value = row;
                VecSetValues(vector, 1, &row, &value, INSERT_VALUES);
            }
        }

        MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(vector);
        VecAssemblyEnd(vector);

        OutputFileHandler handler("DumpPetscObjects");
        std::string output_dir = handler.GetOutputDirectoryFullPath();

        PetscTools::DumpPetscObject(matrix, output_dir + "ten_times_ten.mat");
        PetscTools::DumpPetscObject(vector, output_dir + "ten_times_ten.vec");

        MatDestroy(matrix);
        VecDestroy(vector);

        Mat matrix_read;
        Vec vector_read;

        PetscTools::ReadPetscObject(matrix_read, output_dir + "ten_times_ten.mat");
        PetscTools::ReadPetscObject(vector_read, output_dir + "ten_times_ten.vec");

        double *p_vector_read;
        VecGetArray(vector_read, &p_vector_read);

        for (PetscInt row=0; row<10; row++)
        {
            if (lo<=row && row<hi)
            {
                for (PetscInt col=0; col<10; col++)
                {
                    double value;
                    MatGetValues(matrix_read, 1, &row, 1, &col, &value);
                    TS_ASSERT_EQUALS(value, (double) 10*row+col+1);
                }

            unsigned local_index = row-lo;
            TS_ASSERT_EQUALS(p_vector_read[local_index], (double)row);
            }
        }

        VecRestoreArray(vector_read, &p_vector_read);

        MatDestroy(matrix_read);
        VecDestroy(vector_read);

    }
};
#endif /*TESTPETSCTOOLS_HPP_*/
