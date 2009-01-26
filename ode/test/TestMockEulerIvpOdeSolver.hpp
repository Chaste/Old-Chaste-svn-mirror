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
