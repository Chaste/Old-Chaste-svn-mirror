#ifndef _TESTMOCKEULERIVPODESOLVER_HPP_
#define _TESTMOCKEULERIVPODESOLVER_HPP_
#include <cxxtest/TestSuite.h>

#include <vector>

#include "AbstractIvpOdeSolver.hpp"
#include "MockEulerIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "Ode1.hpp"

class TestMockEulerIvpOdeSolver: public CxxTest::TestSuite
{
    public:
    
    void testMockEulerSolver()
    {
        Ode1 ode_system;
        
        // Initialising the instance of our solver class    
        MockEulerIvpOdeSolver euler_solver;
        // Initialising the instance of our solution class
        OdeSolution solutions;
        
        // Solving the ode problem and writing to solution
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.001, 2.0);
        
        int last = solutions.GetNumberOfTimeSteps();        
        // Test to see if this worked       
        double testvalue = solutions.rGetSolutions()[last][0];
        
        TS_ASSERT_DELTA(testvalue,2.0,0.01);
        
        TS_ASSERT_EQUALS(euler_solver.GetCallCount(), 1);
        
        ode_system.SetInitialConditionsComponent(0,0.0);
        
        state_variables = ode_system.GetInitialConditions();
        solutions = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.001, 2.0);
                                               
        last = solutions.GetNumberOfTimeSteps();        
        // Test to see if this worked       
        testvalue = solutions.rGetSolutions()[last][0];
        
        TS_ASSERT_DELTA(testvalue,2.0,0.01);
        
        TS_ASSERT_EQUALS(euler_solver.GetCallCount(), 2);
        
    }

};
#endif //_TESTMOCKEULERIVPODESOLVER_HPP_
