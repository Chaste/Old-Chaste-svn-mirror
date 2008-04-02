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

    void TestMockEulerSolver()
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
        
        TS_ASSERT_EQUALS(euler_solver.GetCallCount(), 1U);
        
        ode_system.SetInitialConditionsComponent(0,0.0);
        
        state_variables = ode_system.GetInitialConditions();
        solutions = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.001, 2.0);
        
        last = solutions.GetNumberOfTimeSteps();
        // Test to see if this worked
        testvalue = solutions.rGetSolutions()[last][0];
        
        TS_ASSERT_DELTA(testvalue,2.0,0.01);
        
        TS_ASSERT_EQUALS(euler_solver.GetCallCount(), 2U);
        
    }
    
};
#endif //_TESTMOCKEULERIVPODESOLVER_HPP_
