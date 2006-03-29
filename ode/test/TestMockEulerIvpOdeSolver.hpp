#ifndef _TESTMOCKEULERIVPODESOLVER_HPP_
#define _TESTMOCKEULERIVPODESOLVER_HPP_
#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>

#include "AbstractIvpOdeSolver.hpp"
#include "MockEulerIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "TestOde1.hpp"

class TestMockEulerIvpOdeSolver: public CxxTest::TestSuite
{
    public:
    
    void testMockEulerSolver()
    {
        TestOde1 ode_system;
        
        // Initialising the instance of our solver class    
        MockEulerIvpOdeSolver euler_solver;
        // Initialising the instance of our solution class
        OdeSolution solutions;
        
        // Solving the ode problem and writing to solution      
        solutions = euler_solver.Solve(&ode_system, 0.0, 2.0, 0.001, ode_system.GetInitialConditions());
        
        int last = solutions.GetNumberOfTimeSteps();        
        // Test to see if this worked       
        double testvalue = solutions.mSolutions[last][0];
        
        TS_ASSERT_DELTA(testvalue,2.0,0.01);
        
        TS_ASSERT_EQUALS(euler_solver.GetCallCount(), 1);
        
        ode_system.SetInitialConditionsComponent(0,0.0);
        
        solutions = euler_solver.Solve(&ode_system, 0.0, 2.0, 0.001, 
                                       ode_system.GetInitialConditions());
                                               
        last = solutions.GetNumberOfTimeSteps();        
        // Test to see if this worked       
        testvalue = solutions.mSolutions[last][0];
        
        TS_ASSERT_DELTA(testvalue,2.0,0.01);
        
        TS_ASSERT_EQUALS(euler_solver.GetCallCount(), 2);
        
    }

};
#endif //_TESTMOCKEULERIVPODESOLVER_HPP_
