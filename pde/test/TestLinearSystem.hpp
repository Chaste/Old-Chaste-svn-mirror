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
        for (int row=0; row<3; row++)
        {
      	    for(int col=0; col<3; col++)
        	{
        		ls.SetMatrixElement(row, col, (double) row*3+col+1);
       	    }
        }
        ls.AssembleFinalMatrix();
        
        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);
        
        SimpleLinearSolver solver;
        Vec solution_vector;
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve(&solver));
        
        PetscScalar *solution_elements;
        VecGetArray(solution_vector, &solution_elements);
        TS_ASSERT_DELTA(solution_elements[0], 1.0, 0.000001);
        TS_ASSERT_DELTA(solution_elements[1], 2.0, 0.000001);
        TS_ASSERT_DELTA(solution_elements[2], 3.0, 0.000001);
        
        VecDestroy(solution_vector);
        
    }

    void TestLinearSystem2( void )
    {
        LinearSystem ls(2);
        ls.SetMatrixRow(0, 1.0);
        ls.SetMatrixRow(1, 3.0);
        ls.AssembleIntermediateMatrix();
        
        ls.AddToMatrixElement(0, 1, 1.0);
        ls.AddToMatrixElement(1, 1, 1.0);
        ls.AssembleFinalMatrix();

        ls.AddToRhsVectorElement(0, 3.0);
        ls.AddToRhsVectorElement(1, 7.0);
        
        SimpleLinearSolver solver;
        Vec solution_vector;
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve(&solver));
        
        PetscScalar *solution_elements;
        VecGetArray(solution_vector, &solution_elements);
        TS_ASSERT_DELTA(solution_elements[0], 1.0, 0.000001);
        TS_ASSERT_DELTA(solution_elements[1], 1.0, 0.000001);
        
        VecDestroy(solution_vector);
        
    }
};
#endif //_TESTLINEARSYSTEM_HPP_
