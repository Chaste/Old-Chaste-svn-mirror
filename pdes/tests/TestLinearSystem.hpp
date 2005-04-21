#ifndef _TESTLINEARSYSTEM_HPP_
#define _TESTLINEARSYSTEM_HPP_


#include <cxxtest/TestSuite.h>
#include "LinearSystem.hpp"
#include "SimpleLinearSolver.hpp"
          
//class MyTestSuite : public CxxTest::TestSuite 
//{
//public:
//    void testLinearSystem( void )
//    {
//        LinearSystem ls(3);
//        for (int row=0; row<3; row++)
//        {
//      	    for(int col=0; col<3; col++)
//        	{
//        		ls.SetMatrixElement(row, col, (double) row*3+col+1);
//       	    }
//        }
//        ls.AssembleFinalMatrix();
//        
//        ls.SetRhsVectorElement(0, 14.0);
//        ls.SetRhsVectorElement(1, 32.0);
//        ls.SetRhsVectorElement(2, 50.0);
//        
//        SimpleLinearSolver solver;
//        Vec solution_vector = ls.Solve(&solver);
//        
//        PetscScalar *solution_elements;
//        VecGetArray(solution_vector, &solution_elements);
//        TS_ASSERT_DELTA(solution_elements[0], 1.0, 0.000001);
//        TS_ASSERT_DELTA(solution_elements[1], 2.0, 0.000001);
//        TS_ASSERT_DELTA(solution_elements[2], 3.0, 0.000001);
//        
//    }
//};

#endif //_TESTLINEARSYSTEM_HPP_
