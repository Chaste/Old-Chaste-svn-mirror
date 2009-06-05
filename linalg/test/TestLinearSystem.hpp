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

class TestLinearSystem : public CxxTest::TestSuite
{
public:

    void TestLinearSystem1( void )
    {
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

        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);

        // for coverage
        ls.DisplayMatrix();
        ls.DisplayRhs();

        Vec solution_vector;
        solution_vector = ls.Solve();

        int lo,hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        PetscScalar *p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], global_index+1.0, 1e-8);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);
        VecDestroy(solution_vector);

        //SetRelativeTolerance
        ls.SetRelativeTolerance(1e-2);
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve());
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], global_index+1.0, 2e-1);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);
        VecDestroy(solution_vector);

        //SetAbsoluteTolerance
        ls.SetAbsoluteTolerance(1e-5);
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve());
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], global_index+1.0, 1e-8);
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

 
        //Note: these methods are collective.  All processes MUST do them together.
        for (unsigned i=0; i<2; i++)
        {
            ls.ZeroMatrixRow(i);
        }
        
        if (lo <=0 && 0<hi)
        {
            TS_ASSERT_EQUALS(ls.GetMatrixElement(0, 1), 0.0);
        }
        if (lo <=1 && 1<hi)
        {
            TS_ASSERT_EQUALS(ls.GetMatrixElement(1, 1), 0.0);
        }

        if (lo <=2 && 2<hi)
        {
            TS_ASSERT_EQUALS(ls.GetMatrixElement(2, 1), 125.0);
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

    void TestCreateFromVector(void)
    {
        const int SIZE = 5;
        Vec test_vec;
        VecCreate(PETSC_COMM_WORLD, &test_vec);
        VecSetSizes(test_vec, PETSC_DECIDE, SIZE);
        VecSetFromOptions(test_vec);
        LinearSystem ls(test_vec);

        // Check ownership ranges match
        int lo1, hi1, lo2, hi2;
        VecGetOwnershipRange(test_vec, &lo1, &hi1);
        ls.GetOwnershipRange(lo2, hi2);
        TS_ASSERT_EQUALS(lo1, lo2);
        TS_ASSERT_EQUALS(hi1, hi2);

        VecDestroy(test_vec);
    }

    void TestLinearSystem2( void )
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
        PetscScalar *p_solution_elements_array;
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
    void TestNullBasis(void)
    {
        const int SIZE = 5;
        Vec test_vec;
        VecCreate(PETSC_COMM_WORLD, &test_vec);
        VecSetSizes(test_vec, PETSC_DECIDE, SIZE);
        VecSetFromOptions(test_vec);
        LinearSystem ls(test_vec);
        ls.SetNullBasis(&test_vec, 1);

        VecDestroy(test_vec);
    }

    // Test the 3rd constructor
    void TestCreateWithPetscObjects()
    {
        // First try giving it just a vector
        unsigned size = 5u;
        DistributedVector::SetProblemSize(size);
        Vec test_vec = DistributedVector::CreateVec();

        DistributedVector dist_vec(test_vec);
        double test_val = -1.0;
        if (dist_vec.Begin() != dist_vec.End())
        {
            dist_vec[dist_vec.Begin()] = test_val;
        }
        dist_vec.Restore();

        LinearSystem lsv(test_vec, NULL);
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
        PetscTools::SetupMat(m, size, size);

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

    void TestLinearSystem1WithIntialGuess( void )
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
        Vec good_guess;
        VecCreate(PETSC_COMM_WORLD, &good_guess);
        VecSetSizes(good_guess,PETSC_DECIDE,3);
        VecSetFromOptions(good_guess);
        VecSetValue(good_guess, 0, 1.0, INSERT_VALUES);
        VecSetValue(good_guess, 1, 2.0, INSERT_VALUES);
        VecSetValue(good_guess, 2, 3.0, INSERT_VALUES);


        Vec solution_vector;
        solution_vector = ls.Solve(good_guess);
        int lo,hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        PetscScalar *p_solution_elements_array;
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
#if (PETSC_VERSION_MINOR == 2) //Old API
        VecSet(&too_big, bad_guess);
#else
        VecSet(bad_guess, too_big);
#endif
        TS_ASSERT_THROWS_ANYTHING(solution_vector = ls.Solve(bad_guess));

        VecDestroy(solution_vector);
        VecDestroy(good_guess);
        VecDestroy(bad_guess);

    }

    void TestAddMultipleValues( void )
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

        ls.SetMatrixIsSymmetric();

        // Enter symmetric data
        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, fabs(row-col));
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
        double *p_solution_elements_array;
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
        double *p_solution_elements_array, *p_solution_elements_array2;
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
        // others should be their petsc defaults
        TS_ASSERT_EQUALS(atol, 1e-50);
        TS_ASSERT_EQUALS(dtol, 10000.0);
        TS_ASSERT_EQUALS(maxits, 10000);

        KSPType solver;
        PCType pc;
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


    void TestSaveAndLoad()
    {
         //Archive                           
        OutputFileHandler handler("Archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "ls.arch";       
        // SAVE
        {
            LinearSystem ls = LinearSystem(3);
    
            ls.SetMatrixIsSymmetric();
    
            // Enter symmetric data
            for (int row=0; row<3; row++)
            {
                for (int col=0; col<3; col++)
                {
                    ls.SetMatrixElement(row, col, fabs(row-col));
                }
            }
            ls.AssembleFinalLinearSystem();
    
            // arbitrary
            ls.SetRhsVectorElement(0, 14.0);
            ls.SetRhsVectorElement(1, 32.0);
            ls.SetRhsVectorElement(2, 50.0);
                
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            LinearSystem* const p_linear_system = &ls;
            output_arch << p_linear_system;  
            
            TS_ASSERT_EQUALS(p_linear_system->GetSize(), 3u);    
            int lo, hi;
            VecGetOwnershipRange(p_linear_system->GetRhsVector(), &lo, &hi);
            std::vector<double> answer;
            answer.push_back(14.0);
            answer.push_back(32.0);
            answer.push_back(50.0);
            
            for ( int i = lo; i < hi; i++ )
            {       
                TS_ASSERT_DELTA(p_linear_system->GetRhsVectorElement(i), answer[i], 1e-9);
            }

            for (int row=lo; row<hi; row++)
            {
                for (int col=0; col<3; col++)
                {
                    TS_ASSERT_DELTA(p_linear_system->GetMatrixElement(row, col), fabs(row-col), 1e-9);
                }
            }
            
        }
        // LOAD
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs); 
            
            LinearSystem* p_linear_system;
            input_arch >> p_linear_system;
            
            TS_ASSERT_EQUALS(p_linear_system->GetSize(), 3u);    
            
            int lo, hi;
            VecGetOwnershipRange(p_linear_system->GetRhsVector(), &lo, &hi);
            
            std::vector<double> answer;
            answer.push_back(14.0);
            answer.push_back(32.0);
            answer.push_back(50.0);
            
            for ( int i = lo; i < hi; i++ )
            {
                TS_ASSERT_DELTA(p_linear_system->GetRhsVectorElement(i), answer[i], 1e-9);
            }
            
            for (int row=lo; row<hi; row++)
            {
                for (int col=0; col<3; col++)
                {
                    TS_ASSERT_DELTA(p_linear_system->GetMatrixElement(row, col), fabs(row-col), 1e-9);
                }
            }
            
            delete p_linear_system ;
        }
        
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

        KSPType solver;
        PCType pc;
        PC prec;
        KSPGetType(ls.mKspSolver, &solver);
        KSPGetPC(ls.mKspSolver, &prec);
        PCGetType(prec, &pc);

        TS_ASSERT( strcmp(solver,"gmres")==0 );
        TS_ASSERT( strcmp(pc,"jacobi")==0 );
    }
    // the above test should be last in the suite
};
#endif //_TESTLINEARSYSTEM_HPP_
