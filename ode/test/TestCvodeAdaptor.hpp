/*

Copyright (C) University of Oxford, 2008

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

#ifndef _TESTCVODEADAPTOR_HPP_
#define _TESTCVODEADAPTOR_HPP_

#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>
#include <cmath>

#include "CvodeAdaptor.hpp"
#include "AbstractOdeSystem.hpp"
#include "Ode1.hpp"
#include "OdeFirstOrder.hpp"
#include "OdeSecondOrder.hpp"
#include "OdeSecondOrderWithEvents.hpp"


class TestCvodeAdaptor: public CxxTest::TestSuite
{
private:
#ifdef CHASTE_CVODE
    void HelperTestOde1(double startTime, double endTime, double samplingTime)
    {
        // Initialise
        Ode1 ode_system;
        OdeSolution solutions;
        CvodeAdaptor solver;
        
        // Solving the ode problem.
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = solver.Solve(&ode_system, state_variables,
                                 startTime, endTime, endTime-startTime, samplingTime);
        
        int num_timesteps = solutions.GetNumberOfTimeSteps();
        
        // the number of timesteps should be (just about) equal to
        // sim_time/sampling_time
        TS_ASSERT_DELTA(num_timesteps, (endTime-startTime)/samplingTime, 1);
        // also check the size of the data is correct
        TS_ASSERT_EQUALS(solutions.rGetSolutions().size(),
                         (unsigned) (num_timesteps+1));
        
        int last = num_timesteps;
        
        // Test the solution is correct
        double testvalue = solutions.rGetSolutions()[last][0];
        
        // exact solution of Ode1 is y=t-t0
        TS_ASSERT_DELTA(testvalue, endTime-startTime, 0.01);
        
        // Test second version of Solve
        ode_system.SetStateVariables(ode_system.GetInitialConditions());
        state_variables = ode_system.rGetStateVariables();
        solver.Solve(&ode_system, state_variables,
                     startTime, endTime, endTime-startTime);
        TS_ASSERT_DELTA(state_variables[0], endTime-startTime, 0.01);
        
        // no stopping event was specified in the ODE, so check the
        // solver correctly states it didn't stop due to a
        // stopping event.
        TS_ASSERT_EQUALS(solver.StoppingEventOccurred(), false);
    }
#endif // CHASTE_CVODE
    
public:
    void TestBasics() throw (Exception)
    {
#ifdef CHASTE_CVODE
        CvodeAdaptor solver;
        solver.SetMaxSteps(1000);
        TS_ASSERT_EQUALS(solver.GetMaxSteps(), 1000);
        
        TS_ASSERT_DELTA(solver.GetRelativeTolerance(), 1e-4, 1e-12);
        TS_ASSERT_DELTA(solver.GetAbsoluteTolerance(), 1e-6, 1e-12);
        
        solver.SetTolerances(1e-5, 1e-5);
        TS_ASSERT_DELTA(solver.GetRelativeTolerance(), 1e-5, 1e-12);
        TS_ASSERT_DELTA(solver.GetAbsoluteTolerance(), 1e-5, 1e-12);
#endif // CHASTE_CVODE
    }

    void TestOnOde1() throw (Exception)
    {
#ifdef CHASTE_CVODE
    	HelperTestOde1(0.0, 2.0, 0.001);
        HelperTestOde1(1.0, 2.0, 0.01);
        HelperTestOde1(-1.0, 2.0, 2);
        HelperTestOde1(0.0, 0.4, 0.34);
#endif // CHASTE_CVODE
    }
    
    void TestGlobalError() throw (Exception)
    {
#ifdef CHASTE_CVODE
        OdeFirstOrder ode_system;
        
        double h_value=0.01;
        
        CvodeAdaptor solver;
        solver.SetMaxSteps(1000);
        OdeSolution solutions;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = solver.Solve(&ode_system, state_variables,
                                 0.0, 2.0, h_value, 0.1);
        int last = solutions.GetNumberOfTimeSteps();
        double testvalue = solutions.rGetSolutions()[last][0];
        
        // The tests
        double exact_solution = exp(2);

        // TODO: Work out what the global error should be bounded by
        double global_error = 1e-3;

        TS_ASSERT_DELTA(testvalue, exact_solution, global_error);
#endif // CHASTE_CVODE
    }
    
    
    void TestGlobalErrorSystemOf2Equations() throw (Exception)
    {
#ifdef CHASTE_CVODE
        OdeSecondOrder ode_system;
        
        double h_value=0.01;
        
        CvodeAdaptor solver;
        OdeSolution solutions;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = solver.Solve(&ode_system, state_variables,
                                 0.0, 2.0, h_value, 0.1);
        int last = solutions.GetNumberOfTimeSteps();
        
        double testvalue[2];
        testvalue[0] = solutions.rGetSolutions()[last][0];
        testvalue[1] = solutions.rGetSolutions()[last][1];
        
        double exact_solution[2];
        exact_solution[0] = sin(2);
        exact_solution[1] = cos(2);
        
        // TODO: Work this out properly
        double global_error = 1e-3;

        TS_ASSERT_DELTA(testvalue[0],exact_solution[0],global_error);
        TS_ASSERT_DELTA(testvalue[1],exact_solution[1],global_error);
#endif // CHASTE_CVODE
    }
    
    void TestWithStoppingEvent() throw (Exception)
    {
#ifdef CHASTE_CVODE
        OdeSecondOrderWithEvents ode_system;
        CvodeAdaptor solver;
        OdeSolution solutions;
        
        solver.CheckForStoppingEvents();
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.1, 0.01);
        int num_timesteps = solutions.GetNumberOfTimeSteps();
        
        // final time should be about pi/2
        double tend = solutions.rGetTimes()[num_timesteps];
        double tstop = solver.GetStoppingTime();
        TS_ASSERT_DELTA( tend, M_PI_2, 0.01);
        TS_ASSERT_EQUALS(tend, tstop);
        // penultimate y0 should be greater than zero
        TS_ASSERT_LESS_THAN( 0, solutions.rGetSolutions()[num_timesteps-1][0]);
        // final y0 should be less than zero
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[num_timesteps][0], 0);
        // solver should correctly state the stopping event occurred
        TS_ASSERT_EQUALS(solver.StoppingEventOccurred(), true);
        
        // Alternative Solve method
        state_variables = ode_system.GetInitialConditions();
        solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.01);
        TS_ASSERT_EQUALS(solver.StoppingEventOccurred(), true);
        TS_ASSERT_DELTA(solver.GetStoppingTime(), M_PI_2, 0.01)
#endif // CHASTE_CVODE
    }
};

#endif //_TESTCVODEADAPTOR_HPP_
