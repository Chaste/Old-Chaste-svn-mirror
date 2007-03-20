#ifndef _TESTLINEARSYSTEM_HPP_
#define _TESTLINEARSYSTEM_HPP_


#include <cxxtest/TestSuite.h>
#include "LinearSystem.hpp"
#include "SimpleLinearSolver.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestLinearSystem : public CxxTest::TestSuite
{
public:

    void TestLinearSystem1( void )
    {
    
        LinearSystem ls(3);
        
        TS_ASSERT_EQUALS(ls.GetSize(),3U);
        
        
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
    }
};
#endif //_TESTLINEARSYSTEM_HPP_
