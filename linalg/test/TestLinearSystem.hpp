#ifndef _TESTLINEARSYSTEM_HPP_
#define _TESTLINEARSYSTEM_HPP_


#include <cxxtest/TestSuite.h>
#include "LinearSystem.hpp"
#include "SimpleLinearSolver.hpp"
#include "DistributedVector.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include <cmath>


class TestLinearSystem : public CxxTest::TestSuite
{
public:

    void TestLinearSystem1( void )
    {
        LinearSystem ls(3);
        
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
        
        SimpleLinearSolver solver(1e-6);
        Vec solution_vector;
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve(&solver));
        
        int lo,hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        PetscScalar *p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);
        
        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], global_index+1.0, 0.000001);
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
        ls.ZeroLinearSystem();
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
        
        if (lo < 5)
        {
            ls.SetMatrixRow(lo, 1.0);
            ls.AssembleFinalLinearSystem();
            TS_ASSERT_EQUALS(ls.GetMatrixElement(lo, 1), 1.0);
            ls.ZeroMatrixRow(lo);
            ls.AssembleFinalLinearSystem();
            TS_ASSERT_EQUALS(ls.GetMatrixElement(lo, 1), 0.0);
        }
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
        
        SimpleLinearSolver solver(1e-6);
        Vec solution_vector;
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve(&solver));
        
        
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
        
        
        SimpleLinearSolver solver(1e-6);
        
        Vec solution_vector;
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve(&solver, good_guess));
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
        TS_ASSERT_THROWS_ANYTHING(solution_vector = ls.Solve(&solver, bad_guess));
        
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
        
        // solving should be fine
        SimpleLinearSolver solver(1e-6);
        Vec solution_vector;
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve(&solver));
        
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
        SimpleLinearSolver solver(1e-6);
        Vec solution_vector;
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve(&solver));
        

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
        SimpleLinearSolver solver2(1e-6);
        solution_vector2 = ls2.Solve(&solver2);
        
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
    
    void TestSolverCannotHaveDifferentMatricesInput() throw (Exception)
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
        SimpleLinearSolver solver(1e-6);
        Vec solution_vector;
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve(&solver));
	VecDestroy(solution_vector);
        
        //Do it all again
        LinearSystem ls2 = LinearSystem(3);
        
        // Enter non-symmetric data
        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls2.SetMatrixElement(row, col, (double)(10+row-col));
            }
        }
        ls2.AssembleFinalLinearSystem();
        
        // arbitrary
        ls2.SetRhsVectorElement(0, 14.0);
        ls2.SetRhsVectorElement(1, 32.0);
        ls2.SetRhsVectorElement(2, 50.0);
        TS_ASSERT_THROWS_ANYTHING(solution_vector = ls2.Solve(&solver));
    }
};
#endif //_TESTLINEARSYSTEM_HPP_
