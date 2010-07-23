/*

Copyright (C) University of Oxford, 2005-2010

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


#ifndef _TESTLINEARSYSTEM_HPP_
#define _TESTLINEARSYSTEM_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <cmath>
#include "LinearSystem.hpp"
#include "DistributedVector.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "DistributedVectorFactory.hpp"

#include <iostream>
#include <cstring>
#include "ReplicatableVector.hpp"

class TestLinearSystem : public CxxTest::TestSuite
{
public:
   void TestLinearSystem1()
    {
        TS_ASSERT_THROWS_THIS(LinearSystem too_big_to_be_dense(20), "You must provide a rowPreallocation argument for a large sparse system");
        
        LinearSystem ls(3);
        ls.SetMatrixIsConstant(true);

        TS_ASSERT_EQUALS(ls.GetSize(), 3U);

        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double) row*3+col+1);
            }
        }
        ls.AssembleFinalLinearSystem();

        ls.SetRhsVectorElement(0, 1400000.0);
        ls.SetRhsVectorElement(1, 3200000.0);
        ls.SetRhsVectorElement(2, 5000000.0);

        // for coverage
        ls.DisplayMatrix();
        ls.DisplayRhs();

        Vec solution_vector;
        solution_vector = ls.Solve();

        int lo,hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        PetscScalar* p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], 100000.0*(global_index+1), 1e-8);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);
        VecDestroy(solution_vector);

        //SetRelativeTolerance
        ls.SetRelativeTolerance(1e-2);
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve());

        KSPConvergedReason reason;
        KSPGetConvergedReason(ls.mKspSolver, &reason);
        TS_ASSERT_EQUALS(reason, KSP_CONVERGED_RTOL);



        VecGetOwnershipRange(solution_vector,&lo,&hi);
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], 100000.0*(global_index+1), 2e4);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);
        VecDestroy(solution_vector);

        //SetAbsoluteTolerance
        ls.SetAbsoluteTolerance(1e-8);
        solution_vector = ls.Solve();
        KSPGetConvergedReason(ls.mKspSolver, &reason);
        TS_ASSERT_EQUALS(reason, KSP_CONVERGED_ATOL);

        //Check that it converged for the right reason

        VecGetOwnershipRange(solution_vector,&lo,&hi);
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], 100000.0*(global_index+1), 1e-8);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);
        VecDestroy(solution_vector);
    }

    void TestZeroingLinearSystem()
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixRow(i, (double)i);
            ls.SetRhsVectorElement(i, (double)i);
        }
        ls.AssembleFinalLinearSystem();

        int lo, hi;
        ls.GetOwnershipRange(lo, hi);
        for (int row=lo; row<hi; row++)
        {
            TS_ASSERT_EQUALS(ls.GetRhsVectorElement(row), row);
            for (int col=0; col<5; col++)
            {
                TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), row);
            }
        }

        ls.ZeroLinearSystem();

        for (int row=lo; row<hi; row++)
        {
            TS_ASSERT_EQUALS(ls.GetRhsVectorElement(row), 0.0);
            for (int col=0; col<5; col++)
            {
                TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), 0.0);
            }
        }


        ls.SetMatrixRow(2, 125.0);

        //Note: this method is collective.  All processes MUST do it together.
        ls.AssembleFinalLinearSystem();


        if (lo <=2 && 2<hi)
        {
            TS_ASSERT_EQUALS(ls.GetMatrixElement(2, 1), 125.0);
        }
    }

    void TestZeroMatrixRowsWithValueOnDiagonal()
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixRow(i, (double)i);
            ls.SetRhsVectorElement(i, (double)i);
        }
        ls.AssembleFinalLinearSystem();

        std::vector<unsigned> rows;
        rows.push_back(2);
        rows.push_back(3);
        rows.push_back(4);

        ls.ZeroMatrixRowsWithValueOnDiagonal(rows, 3.14);

        int lo, hi;
        ls.GetOwnershipRange(lo, hi);
        for(int row=2; row<5; row++)
        {
            if (lo<=row && row<hi)
            {
                for(int i=0; i<5; (i+1==row? i+=2 : i++)) // for i=0,1..,row-1,row+1,..,5
                {
                    TS_ASSERT_EQUALS(ls.GetMatrixElement(row,i), 0.0);
                }
                TS_ASSERT_EQUALS(ls.GetMatrixElement(row,row), 3.14);
            }
        }
    }



    void TestZeroingLinearSystemByColumn()
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixElement(i, i, 3.0);
        }
        ls.SetMatrixElement(0, 1, 4.0);

        ls.AssembleFinalLinearSystem();

        for (unsigned col=0; col<5; col++)
        {
            ls.ZeroMatrixColumn(col);
        }
        ls.AssembleFinalLinearSystem();

        int lo, hi;
        ls.GetOwnershipRange(lo, hi);
        for (int row=lo; row<hi; row++)
        {
            TS_ASSERT_EQUALS(ls.GetRhsVectorElement(row), 0);
            for (int col=0; col<5; col++)
            {
                TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), 0);
            }
        }

        MatInfo info;
        double num_nonzeros;

        MatGetInfo(ls.rGetLhsMatrix(),MAT_GLOBAL_SUM,&info);

        num_nonzeros = info.nz_used;

        TS_ASSERT_EQUALS(int(num_nonzeros),6);
    }

    void TestZeroMatrixRowsAndColumnsWithValueOnDiagonal() throw(Exception)
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixRow(i, (double)i);
        }
        ls.AssembleFinalLinearSystem();

        std::vector<unsigned> rows(3);
        rows[0] = 2;
        rows[1] = 3;
        rows[2] = 4;

        ls.ZeroMatrixRowsAndColumnsWithValueOnDiagonal(rows, 3.1);

        int lo, hi;
        ls.GetOwnershipRange(lo, hi);
        for (int row=lo; row<hi; row++)
        {
            for (int col=0; col<5; col++)
            {
                if( (col>=2) || (row>=2) )
                {
                    // the altered values
                    if(row!=col)
                    {
                        TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), 0);
                    }
                    else
                    {
                        TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), 3.1);
                    }
                }
                else
                {
                    // unaltered values
                    TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), row);
                }
            }
        }
    }

    void TestGetMatrixRowDistributed()
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixRow(i, (double)i);
        }
        ls.AssembleFinalLinearSystem();

        Vec third_row = ls.GetMatrixRowDistributed(3);

        DistributedVectorFactory factory(third_row);

        DistributedVector distributed_third_row = factory.CreateDistributedVector(third_row);
        for (DistributedVector::Iterator index = distributed_third_row.Begin();
             index!= distributed_third_row.End();
             ++index)
        {
            TS_ASSERT_EQUALS(distributed_third_row[index], 3);
        }

        VecDestroy(third_row);
    }


    void TestCreateFromVector()
    {
        const int SIZE = 5;
        Vec test_vec=PetscTools::CreateVec(SIZE);

        LinearSystem ls(test_vec, 5);

        // Check ownership ranges match
        int lo1, hi1, lo2, hi2;
        VecGetOwnershipRange(test_vec, &lo1, &hi1);
        ls.GetOwnershipRange(lo2, hi2);
        TS_ASSERT_EQUALS(lo1, lo2);
        TS_ASSERT_EQUALS(hi1, hi2);

        VecDestroy(test_vec);
    }

    void TestLinearSystem2()
    {
        LinearSystem ls(2);
        ls.SetMatrixRow(0, 1.0);
        ls.SetMatrixRow(1, 3.0);
        ls.AssembleIntermediateLinearSystem();

        ls.AddToMatrixElement(0, 1, 1.0);
        ls.AddToMatrixElement(1, 1, 1.0);
        ls.AssembleFinalLinearSystem();

        ls.AddToRhsVectorElement(0, 3.0);
        ls.AddToRhsVectorElement(1, 7.0);

        Vec solution_vector;
        solution_vector = ls.Solve();


        int lo,hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        PetscScalar* p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<2; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], 1.0, 0.000001);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);

        VecDestroy(solution_vector);
    }

    /**
     * This is a stub test for coverage purposes.
     */
    void TestNullBasis()
    {
 #ifndef NDEBUG //Only do this test in debug mode, since functionality is skipped in optimized code
        unsigned size = 10;

        // Test it throws if one of the vectors in the base is not normal
        {
            std::vector<double> data(size,1.0);
            Vec non_orthonormal = PetscTools::CreateVec(data);

            LinearSystem ls((PetscInt) size);
            TS_ASSERT_THROWS_THIS(ls.SetNullBasis(&non_orthonormal, 1),
                    "One of the vectors in the null space is not normalised");
            VecDestroy(non_orthonormal);
        }

        // Test it throws if the vectors in the base are not orthogonal
        {
            std::vector<double> data(size,0.0);
            data[0] = 1.0;
            Vec one_zeros = PetscTools::CreateVec(data);

            std::vector<double> data2(size,0.0);
            data2[1] = 1.0;
            Vec zero_one_zeros = PetscTools::CreateVec(data2);

            Vec null_basis[] = {one_zeros, zero_one_zeros, one_zeros};

            LinearSystem ls((PetscInt) size);
            TS_ASSERT_THROWS_THIS(ls.SetNullBasis(null_basis, 3),"The null space is not orthogonal.");

            VecDestroy(one_zeros);
            VecDestroy(zero_one_zeros);
        }


        // Test it doesn't throws if a non-orthonormal basis is passed
        {
            std::vector<double> data(size,0.0);
            data[0] = 1.0;
            Vec one_zeros = PetscTools::CreateVec(data);

            std::vector<double> data2(size,0.0);
            data2[1] = 1.0;
            Vec zero_one_zeros = PetscTools::CreateVec(data2);

            Vec null_basis[] = {one_zeros, zero_one_zeros};

            LinearSystem ls((PetscInt) size);
            TS_ASSERT_THROWS_NOTHING(ls.SetNullBasis(null_basis, 2));

            VecDestroy(one_zeros);
            VecDestroy(zero_one_zeros);
        }
#endif
    }
    
   void TestRemoveNullSpace()
    {
        LinearSystem ls(3);
        ls.SetMatrixIsConstant(true);

        TS_ASSERT_EQUALS(ls.GetSize(), 3U);

        for (int row=0; row<3; row++)
        {
            ls.SetMatrixElement(row, row, (double) row+1);
        }
        ls.AssembleFinalLinearSystem();

        ls.SetRhsVectorElement(0, 1.0);
        ls.SetRhsVectorElement(1, 2.0);
        ls.SetRhsVectorElement(2, 3.0);

        std::vector<double> data(3,0.0);
        data[0] = 1.0;
        Vec one_zeros = PetscTools::CreateVec(data);

        Vec null_basis[] = {one_zeros};

        ls.SetNullBasis(null_basis, 1);

        Vec wrong_solution = ls.Solve();
        ReplicatableVector replicated_wrong_solution(wrong_solution);

        // Wrong solution since wrong null space was provided.
        TS_ASSERT_DELTA(replicated_wrong_solution[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(replicated_wrong_solution[1], 1.0, 1e-8);
        TS_ASSERT_DELTA(replicated_wrong_solution[2], 1.0, 1e-8);
        
        // Now remove the null space and we will hopefully get the right solution
        ls.RemoveNullSpace();
        Vec solution = ls.Solve();
        ReplicatableVector replicated_solution(solution);

        TS_ASSERT_DELTA(replicated_solution[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(replicated_solution[1], 1.0, 1e-8);
        TS_ASSERT_DELTA(replicated_solution[2], 1.0, 1e-8);
        
        VecDestroy(one_zeros);
        VecDestroy(wrong_solution);
        VecDestroy(solution);

    }    

    // Test the 3rd constructor
    void TestCreateWithPetscObjects()
    {
        // First try giving it just a vector
        unsigned size = 5u;
        DistributedVectorFactory factory(size);
        Vec test_vec = factory.CreateVec();

        DistributedVector dist_vec = factory.CreateDistributedVector(test_vec);
        double test_val = -1.0;
        if (dist_vec.Begin() != dist_vec.End())
        {
            dist_vec[dist_vec.Begin()] = test_val;
        }
        dist_vec.Restore();

        LinearSystem lsv(test_vec, (Mat) NULL);
        TS_ASSERT_EQUALS(lsv.GetSize(), size);

        if (dist_vec.Begin() != dist_vec.End())
        {
            TS_ASSERT_EQUALS(lsv.GetRhsVectorElement(dist_vec.Begin().Global),
                             test_val);
        }

        // Change the Vec and see if the linear system reflects the change
        double test_val2 = 2.0;
        if (dist_vec.Begin() != dist_vec.End())
        {
            dist_vec[dist_vec.Begin()] = test_val2;
        }
        dist_vec.Restore();
        if (dist_vec.Begin() != dist_vec.End())
        {
            TS_ASSERT_EQUALS(lsv.GetRhsVectorElement(dist_vec.Begin().Global),
                             test_val2);
        }

        // Now try with just a matrix
        Mat m;
        PetscTools::SetupMat(m, size, size, size);

        if (dist_vec.Begin() != dist_vec.End())
        {
            MatSetValue(m, dist_vec.Begin().Global, 0, test_val, INSERT_VALUES);
        }

        LinearSystem lsm(NULL, m);
        TS_ASSERT_EQUALS(lsm.GetSize(), size);
        lsm.AssembleFinalLhsMatrix();

        if (dist_vec.Begin() != dist_vec.End())
        {
            TS_ASSERT_EQUALS(lsm.GetMatrixElement(dist_vec.Begin().Global, 0),
                             test_val);
        }

        // Change the Mat and see if the linear system reflects the change
        if (dist_vec.Begin() != dist_vec.End())
        {
            MatSetValue(m, dist_vec.Begin().Global, 0, test_val2, INSERT_VALUES);
        }
        MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
        if (dist_vec.Begin() != dist_vec.End())
        {
            TS_ASSERT_EQUALS(lsm.GetMatrixElement(dist_vec.Begin().Global, 0),
                             test_val2);
        }

        VecDestroy(test_vec);
        MatDestroy(m);
    }

    void TestLinearSystem1WithIntialGuess()
    {
        LinearSystem ls(3);


        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double) row*3+col+1);
            }
        }
        ls.AssembleFinalLinearSystem();

        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);

        //Set the correct answer for the intial guess
        Vec good_guess=PetscTools::CreateVec(3);
        VecSetValue(good_guess, 0, 1.0, INSERT_VALUES);
        VecSetValue(good_guess, 1, 2.0, INSERT_VALUES);
        VecSetValue(good_guess, 2, 3.0, INSERT_VALUES);


        Vec solution_vector;
        solution_vector = ls.Solve(good_guess);
        int lo,hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        PetscScalar* p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_EQUALS(p_solution_elements_array[local_index], global_index+1.0);
                //Zero tolerance
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);

        //Set the a bad intial guess
        Vec bad_guess;
        VecDuplicate(good_guess, &bad_guess);
        PetscScalar too_big = 1e5;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        VecSet(&too_big, bad_guess);
#else
        VecSet(bad_guess, too_big);
#endif
        TS_ASSERT_THROWS_CONTAINS(solution_vector = ls.Solve(bad_guess),
                "DIVERGED_DTOL in function");
        VecDestroy(solution_vector);
        VecDestroy(good_guess);
        VecDestroy(bad_guess);

    }

    void TestAddMultipleValues()
    {

        LinearSystem syst = LinearSystem(3);

        c_matrix<double, 2, 2> small_matrix;
        c_vector<double, 2> small_vector;

        small_matrix(0,0) = 1;
        small_matrix(0,1) = 2;
        small_matrix(1,0) = 3;
        small_matrix(1,1) = 4;

        small_vector(0) = -1;
        small_vector(1) = -2;

        unsigned large_matrix_indices[2]={0,2};

        syst.AddLhsMultipleValues(large_matrix_indices, small_matrix);
        syst.AddRhsMultipleValues(large_matrix_indices, small_vector);

        syst.AssembleFinalLinearSystem();

        PetscInt lo, hi;
        syst.GetOwnershipRange(lo, hi);

        if (lo <=0 && 0<hi)
        {
            TS_ASSERT_EQUALS(syst.GetMatrixElement(0,0), 1);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(0,1), 0);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(0,2), 2);

            TS_ASSERT_EQUALS(syst.GetRhsVectorElement(0), -1);
        }
        if (lo <=1 && 1<hi)
        {
            TS_ASSERT_EQUALS(syst.GetMatrixElement(1,0), 0);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(1,1), 0);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(1,2), 0);

            TS_ASSERT_EQUALS(syst.GetRhsVectorElement(1), 0);
        }
        if (lo <=2 && 2<hi)
        {
            TS_ASSERT_EQUALS(syst.GetMatrixElement(2,0), 3);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(2,1), 0);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(2,2), 4);

            TS_ASSERT_EQUALS(syst.GetRhsVectorElement(2), -2);
        }
    }


    void TestSymmetricMatrix()
    {
        LinearSystem ls = LinearSystem(3);

        TS_ASSERT(!ls.IsMatrixSymmetric());
        ls.SetMatrixIsSymmetric();
        TS_ASSERT(ls.IsMatrixSymmetric());

        // Enter symmetric data
        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double)(abs(row-col)));
            }
        }
        ls.AssembleFinalLinearSystem();

        // arbitrary
        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);

        Vec solution_vector;
        solution_vector = ls.Solve();

        double expected_solution[3]={25.0,0.0,7.0};
        PetscInt lo, hi;
        ls.GetOwnershipRange(lo, hi);
        double* p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);
        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], expected_solution[global_index], 1e-5);
            }
        }

        VecRestoreArray(solution_vector, &p_solution_elements_array);

        VecDestroy(solution_vector);

        // coverage
        ls.SetMatrixIsSymmetric(false);
        TS_ASSERT(!ls.IsMatrixSymmetric());
    }

    void TestNonSymmetricMatrix()
    {
        LinearSystem ls = LinearSystem(3);

        // Enter non-symmetric data
        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double)(10+row-col));
            }
        }
        ls.AssembleFinalLinearSystem();

        // arbitrary
        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);

        // solving should be fine
        Vec solution_vector;
        solution_vector = ls.Solve();


        LinearSystem ls2 = LinearSystem(3);
        ls2.SetMatrixIsSymmetric();

        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls2.SetMatrixElement(row, col, (double)(10+row-col));
            }
        }
        ls2.AssembleFinalLinearSystem();

        // what happens when we solve?
        Vec solution_vector2;
        solution_vector2 = ls2.Solve();

        //Check answers
        double expected_solution[3]={-68.0,6.0,80.0};
        PetscInt lo, hi;
        ls.GetOwnershipRange(lo, hi);
        double* p_solution_elements_array,* p_solution_elements_array2;
        VecGetArray(solution_vector, &p_solution_elements_array);
        VecGetArray(solution_vector2, &p_solution_elements_array2);
        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], expected_solution[global_index], 1e-5);
                TS_ASSERT_LESS_THAN(2, fabs(p_solution_elements_array2[local_index] - expected_solution[global_index]));
                //Diverges from expected by more than 2
            }
        }

        VecRestoreArray(solution_vector, &p_solution_elements_array);
        VecRestoreArray(solution_vector2, &p_solution_elements_array2);
        VecDestroy(solution_vector2);
        VecDestroy(solution_vector);


    }

    void TestGetSetKSP() throw (Exception)
    {
        /////////////////////////
        // Set relative tolerance before first solve
        /////////////////////////
        LinearSystem ls = LinearSystem(5);
        ls.SetRelativeTolerance(1e-3);
        ls.SetKspType("cg");
        ls.SetPcType("jacobi");
        ls.AssembleFinalLinearSystem();
        Vec solution_vector;
        solution_vector = ls.Solve();
        VecDestroy(solution_vector);
        PetscReal rtol, atol, dtol;
        int maxits;
        KSPGetTolerances(ls.mKspSolver, &rtol, &atol, &dtol, &maxits);
        TS_ASSERT_EQUALS(rtol, 1e-3);
        // others should be their PETSc defaults (unless we've done different)
        TS_ASSERT_EQUALS(atol, 1e-50);
        TS_ASSERT_EQUALS(dtol, 10000.0);
        TS_ASSERT_EQUALS(maxits, 10000);

#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        const KSPType solver;
        const PCType pc;
#else
        KSPType solver;
        PCType pc;
#endif

        PC prec;
        KSPGetType(ls.mKspSolver, &solver);
        KSPGetPC(ls.mKspSolver, &prec);
        PCGetType(prec, &pc);
        TS_ASSERT( strcmp(solver,"cg")==0 );
        TS_ASSERT( strcmp(pc,"jacobi")==0 );
        ls.SetKspType("gmres");
        ls.SetPcType("ilu");
        //Test that we can change the solver type after its first use
        KSPGetType(ls.mKspSolver, &solver);
        PCGetType(prec, &pc);
        TS_ASSERT( strcmp(solver,"gmres")==0 );
        TS_ASSERT( strcmp(pc,"ilu")==0 );


        /////////////////////////////////
        // Set relative tolerance after first solve
        /////////////////////////////////
        ls.SetRelativeTolerance(1e-4);
        KSPGetTolerances(ls.mKspSolver, &rtol, &atol, &dtol, &maxits);
        TS_ASSERT_EQUALS(rtol, 1e-4);
        TS_ASSERT_EQUALS(atol, 1e-50);
        TS_ASSERT_EQUALS(dtol, 10000.0);
        TS_ASSERT_EQUALS(maxits, 10000);

        /////////////////////////////////
        // Set abs tolerance before first solve
        //////////////////////////////////
        LinearSystem ls2 = LinearSystem(5);
        ls2.SetAbsoluteTolerance(1e-3);
        ls2.AssembleFinalLinearSystem();
        Vec solution_vector2;
        solution_vector2 = ls2.Solve();
        VecDestroy(solution_vector2);
        KSPGetTolerances(ls2.mKspSolver, &rtol, &atol, &dtol, &maxits);
        TS_ASSERT_EQUALS(rtol, DBL_EPSILON);
        TS_ASSERT_EQUALS(atol, 1e-3);
        TS_ASSERT_EQUALS(dtol, 10000.0);
        TS_ASSERT_EQUALS(maxits, 10000);

        ///////////////////////////////////
        // Set abs tolerance after first solve
        ////////////////////////////////////
        ls2.SetAbsoluteTolerance(1e-2);
        Vec solution_vector3;
        solution_vector3 = ls2.Solve();
        VecDestroy(solution_vector3);
        KSPGetTolerances(ls2.mKspSolver, &rtol, &atol, &dtol, &maxits);
        TS_ASSERT_EQUALS(rtol, DBL_EPSILON);
        TS_ASSERT_EQUALS(atol, 1e-2);
        TS_ASSERT_EQUALS(dtol, 10000.0);
        TS_ASSERT_EQUALS(maxits, 10000);

    }

    void TestPetscSaveAndLoad()
    {
         //Archive
        OutputFileHandler handler("Archive", false);
        std::string archive_filename_lhs, archive_filename_rhs;
        archive_filename_lhs = handler.GetOutputDirectoryFullPath() + "direct_lhs.mat";
        archive_filename_rhs = handler.GetOutputDirectoryFullPath() + "direct_rhs.vec";

        // Make a linear system
        LinearSystem ls = LinearSystem(3);
        ls.SetMatrixIsSymmetric();
        // Enter symmetric data
        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double)(abs(row-col)));
            }
        }
        ls.AssembleFinalLinearSystem();
        // arbitrary
        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);

        // SAVE
        {
            PetscViewer vec_viewer;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
            PetscViewerFileType type = PETSC_FILE_CREATE;
#else
            PetscFileMode type = FILE_MODE_WRITE;
#endif
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_rhs.c_str(), type, &vec_viewer);

            VecView(ls.GetRhsVector(), vec_viewer);
            PetscViewerDestroy(vec_viewer);

            PetscViewer mat_viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_lhs.c_str(), type, &mat_viewer);

            MatView(ls.GetLhsMatrix(), mat_viewer);
            PetscViewerDestroy(mat_viewer);
        }
        // LOAD
        {

            PetscViewer vec_viewer;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
            PetscViewerFileType type = PETSC_FILE_RDONLY;
#else
            PetscFileMode type = FILE_MODE_READ;
#endif
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_rhs.c_str(), type, &vec_viewer);
            Vec new_vec;
            VecLoad(vec_viewer, PETSC_NULL, &new_vec);
            PetscViewerDestroy(vec_viewer);

            int lo, hi;
            VecGetOwnershipRange(new_vec, &lo, &hi);
            std::vector<double> answer;
            answer.push_back(14.0);
            answer.push_back(32.0);
            answer.push_back(50.0);

            double* p_vec_values;
            VecGetArray(new_vec, &p_vec_values);

            for ( int i = lo; i < hi; i++ )
            {
                TS_ASSERT_DELTA(p_vec_values[i-lo], answer[i], 1e-9);
            }

            VecDestroy(new_vec);

            PetscViewer mat_viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_lhs.c_str(), type, &mat_viewer);
            Mat new_mat;
            MatLoad(mat_viewer, PETSC_NULL, &new_mat);
            PetscViewerDestroy(mat_viewer);

            for (int row=lo; row<hi; row++)
            {
                // Get a whole row out of the matrix and check it
                PetscInt row_as_array[1];
                row_as_array[0] = row;
                PetscInt col_as_array[3];
                for (int col=0; col<3; col++)
                {
                   col_as_array[col] = col;
                }
                double ret_array[3];
                MatGetValues(new_mat, 1, row_as_array, 3, col_as_array, ret_array);

                for (int col=0; col<3; col++)
                {
                    TS_ASSERT_DELTA(ret_array[col], (double)(abs(row-col)), 1e-9);
                }
            }

            MatDestroy(new_mat);
        }

    }

    void TestSaveAndLoadLinearSystem()
    {
         //Archive
        OutputFileHandler handler("Archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "linear_system.arch";

        int lo, hi;
        unsigned size=5;
        std::vector<double> rhs_values;
        rhs_values.push_back(14.0);
        rhs_values.push_back(32.0);
        rhs_values.push_back(50.0);
        rhs_values.push_back(50.0);
        rhs_values.push_back(50.0);
        TS_ASSERT_EQUALS(rhs_values.size(), size);
        // SAVE
        {
            LinearSystem ls = LinearSystem(size);
            Mat temp_mat=ls.GetLhsMatrix();
            PetscTruth symm_set, is_symmetric;
            is_symmetric = PETSC_FALSE;
            MatIsSymmetricKnown(temp_mat, &symm_set, &is_symmetric);
            TS_ASSERT_EQUALS(symm_set, PETSC_FALSE);
            TS_ASSERT_EQUALS(is_symmetric, PETSC_FALSE);

            ls.SetMatrixIsSymmetric();

            MatIsSymmetricKnown(temp_mat, &symm_set, &is_symmetric);
            TS_ASSERT_EQUALS(symm_set, PETSC_TRUE);
            TS_ASSERT_EQUALS(is_symmetric, PETSC_TRUE);

            // Enter symmetric data
            for (unsigned row=0; row<size; row++)
            {
                for (unsigned col=0; col<size; col++)
                {
                    ls.SetMatrixElement(row, col, (row == col)?(row+1.0):0.0);
                }
            }
            ls.AssembleFinalLinearSystem();

            for (unsigned i=0; i<size; i++)
            {
                ls.SetRhsVectorElement(i, rhs_values[i]);
            }
            ls.SetKspType("cg");
            ls.SetPcType("none");

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            LinearSystem* const p_linear_system = &ls;
            output_arch << p_linear_system;

            TS_ASSERT_EQUALS(p_linear_system->GetSize(), size);
            VecGetOwnershipRange(p_linear_system->GetRhsVector(), &lo, &hi);

            for ( int i = lo; i < hi; i++ )
            {
                TS_ASSERT_DELTA(p_linear_system->GetRhsVectorElement(i), rhs_values[i], 1e-9);
            }

            for (unsigned row=(unsigned)lo; row<(unsigned)hi; row++)
            {
                for (unsigned col=0; col<size; col++)
                {
                    TS_ASSERT_DELTA(p_linear_system->GetMatrixElement(row, col), (row == col)?(row+1.0):0.0, 1e-9);
                }
            }

        }
        // LOAD
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            //LinearSystem linear_system(3);
            LinearSystem* p_linear_system;//=&linear_system;
            input_arch >> p_linear_system;

            //Check that structural symmetry is preserved
            PetscTruth symm_set, is_symmetric;
            is_symmetric=PETSC_FALSE;
            MatIsSymmetricKnown(p_linear_system->GetLhsMatrix(), &symm_set, &is_symmetric);
            TS_ASSERT_EQUALS(symm_set, PETSC_TRUE);
            TS_ASSERT_EQUALS(is_symmetric, PETSC_TRUE);


            TS_ASSERT_EQUALS(p_linear_system->GetSize(), size);

            int saved_lo, saved_hi;
            VecGetOwnershipRange(p_linear_system->GetRhsVector(), &saved_lo, &saved_hi);

            TS_ASSERT_EQUALS(hi, saved_hi);
            TS_ASSERT_EQUALS(lo, saved_lo);

            for ( int i = lo; i < hi; i++ )
            {
                TS_ASSERT_DELTA(p_linear_system->GetRhsVectorElement(i), rhs_values[i], 1e-9);
            }

            for (unsigned row=(unsigned)lo; row<(unsigned)hi; row++)
            {
                for (unsigned col=0; col<size; col++)
                {
                    TS_ASSERT_DELTA(p_linear_system->GetMatrixElement(row, col), (row == col)?(row+1.0):0.0, 1e-9);
                }
            }

            //Check archiving of KSP/PC types
            Vec solution_vector3;
            solution_vector3 = p_linear_system->Solve();
            VecDestroy(solution_vector3);
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
            const KSPType solver;
            const PCType pc;
#else
            KSPType solver;
            PCType pc;
#endif
            PC prec;
            KSPGetType(p_linear_system->mKspSolver, &solver);
            KSPGetPC(p_linear_system->mKspSolver, &prec);
            PCGetType(prec, &pc);

            TS_ASSERT( strcmp(solver, "cg")==0 );
            TS_ASSERT( strcmp(pc, "none")==0 );
            delete p_linear_system;
        }

    }

    void TestConsecutiveSolvesDifferentPreconditioner()
    {
            Mat system_matrix;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat");

            Vec system_rhs;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec");

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("cg");

            ls.SetPcType("bjacobi");
            Vec solution = ls.Solve(/*no guess provided*/);
            unsigned block_jacobi_its = ls.GetNumIterations();
            VecDestroy(solution);

            ls.SetPcType("blockdiagonal");
            solution = ls.Solve(/*no guess provided*/);
            unsigned ldu_its = ls.GetNumIterations();
            VecDestroy(solution);

            ls.SetPcType("bjacobi");
            solution = ls.Solve(/*no guess provided*/);
            unsigned second_block_jacobi_its = ls.GetNumIterations();
            VecDestroy(solution);

            TS_ASSERT_DIFFERS(block_jacobi_its, ldu_its)
            TS_ASSERT_DIFFERS(ldu_its, second_block_jacobi_its)

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
    }


    // this test should be the last in the suite
    void TestSetFromOptions()
    {
        LinearSystem ls = LinearSystem(5);
        PetscOptionsSetValue("-ksp_type", "gmres");
        PetscOptionsSetValue("-pc_type", "jacobi");
        ls.AssembleFinalLinearSystem();

        ls.SetKspType("cg"); //Not really -- see above
        ls.SetPcType("ilu"); //Not really -- see above
        Vec solution_vector3;
        solution_vector3 = ls.Solve();
        VecDestroy(solution_vector3);

#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        const KSPType solver;
        const PCType pc;
#else
        KSPType solver;
        PCType pc;
#endif
        PC prec;
        KSPGetType(ls.mKspSolver, &solver);
        KSPGetPC(ls.mKspSolver, &prec);
        PCGetType(prec, &pc);

        TS_ASSERT( strcmp(solver,"gmres")==0 );
        TS_ASSERT( strcmp(pc,"jacobi")==0 );
    }
    // the above test should be last in the suite

//    void TestSingularSolves() throw(Exception)
//    {
//        LinearSystem ls(2);
//
//        ls.SetMatrixElement(0, 0, 2);
//        ls.SetMatrixElement(0, 1, 2);
//        ls.SetMatrixElement(1, 0, 2);
//        ls.SetMatrixElement(1, 1, 2);
//
//        ls.SetRhsVectorElement(0, 100.0);
//        ls.SetRhsVectorElement(1, 100.0);
//
//        ls.AssembleFinalLinearSystem();
//
//        Vec x = ls.Solve();
//        ReplicatableVector xx(x);
//
//        std::cout << xx[0] << " " << xx[1] << "\n"; //solves fine without null space?
//    }
};
#endif //_TESTLINEARSYSTEM_HPP_
