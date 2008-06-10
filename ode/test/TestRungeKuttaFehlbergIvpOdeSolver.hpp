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


#ifndef TESTRUNGEKUTTAFEHLBERGIVPODESOLVER_HPP_
#define TESTRUNGEKUTTAFEHLBERGIVPODESOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "RkfTestOde.hpp"
#include "OdeThirdOrder.hpp"
#include "OdeThirdOrderWithEvents.hpp"
#include "Ode4.hpp"
#include "Ode5.hpp"
#include "Ode5Jacobian.hpp"
#include "VanDerPolOde.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "OutputFileHandler.hpp"

class TestRungeKuttaFehlbergIvpOdeSolver: public CxxTest::TestSuite
{
public:
    void TestAdjustStepSize() throw(Exception)
    {
         RungeKuttaFehlbergIvpOdeSolver rkf_solver;

         double current_step_size = 0.1;
         double error = 1e-04;
         double tolerance = 1e-06;
         double max_step_size = 1.0;
         double min_step_size = 0.002;

         // Use default step size calculation
         double step_size_should_be = current_step_size*pow(tolerance/(2*error),0.25);
         rkf_solver.AdjustStepSize(current_step_size, error,
                                   tolerance, max_step_size, min_step_size);

         TS_ASSERT_DELTA(current_step_size, step_size_should_be, 1e-7);

         // Make step size 1/10th
         error = 1e-1;
         step_size_should_be = 0.1*current_step_size;
         rkf_solver.AdjustStepSize(current_step_size, error,
                                   tolerance, max_step_size, min_step_size);
         TS_ASSERT_DELTA(current_step_size, step_size_should_be, 1e-9);

         // Make step size 4 times bigger
         error = 1e-19;
         step_size_should_be = 4*current_step_size;
         rkf_solver.AdjustStepSize(current_step_size, error,
                                   tolerance, max_step_size, min_step_size);
         TS_ASSERT_DELTA(current_step_size, step_size_should_be, 1e-9);

         // Ramp up step size up to maximum
         error = 1e-19;
         step_size_should_be = max_step_size;
         for (unsigned i=0; i<4; i++)
         {
            // Multiplies step size by 4 each time
            rkf_solver.AdjustStepSize(current_step_size, error,
                                      tolerance, max_step_size, min_step_size);
         }

         TS_ASSERT_DELTA(current_step_size, step_size_should_be, 1e-9);

         // try to make step size too small
         error = 1e10;
         for (unsigned i=0; i<2; i++)
         {
            // Divides step size by 10 each time (still within acceptable bounds)
            rkf_solver.AdjustStepSize(current_step_size, error,
                                      tolerance, max_step_size, min_step_size);
         }
         // Divides step size by 10 and causes it to fall over.
         TS_ASSERT_THROWS_ANYTHING(rkf_solver.AdjustStepSize(current_step_size, error,
                                   tolerance, max_step_size, min_step_size) );
    }

    void TestCalculateNextYValue() throw(Exception)
    {
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;

        RkfTestOde ode;

        double time_step = 0.25;
        double time = 0.0;
        std::vector<double> current_y_values(1);
        std::vector<double> next_y_values(1);

        current_y_values[0] = 0.5;
        next_y_values = current_y_values;

        unsigned number_of_variables = ode.GetNumberOfStateVariables();

        /*
         * These commands are all in the InternalSolve() method, which is
         * the only thing that ever calls CalculateYValues.
         */
        rkf_solver.mError.reserve(number_of_variables);
        rkf_solver.mk1.reserve(number_of_variables);
        rkf_solver.mk2.reserve(number_of_variables);
        rkf_solver.mk3.reserve(number_of_variables);
        rkf_solver.mk4.reserve(number_of_variables);
        rkf_solver.mk5.reserve(number_of_variables);
        rkf_solver.mk6.reserve(number_of_variables);
        rkf_solver.myk2.reserve(number_of_variables);
        rkf_solver.myk3.reserve(number_of_variables);
        rkf_solver.myk4.reserve(number_of_variables);
        rkf_solver.myk5.reserve(number_of_variables);
        rkf_solver.myk6.reserve(number_of_variables);


        rkf_solver.CalculateNextYValue(&ode, time_step, time,
                                        current_y_values, next_y_values);

        TS_ASSERT_DELTA(current_y_values[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(next_y_values[0], 9.204873e-01, 1e-5);
        TS_ASSERT_DELTA(rkf_solver.mError[0], 6.2e-6, 1e-7);
    }

    void TestRKFehlbergWithExampleFromBook() throw(Exception)
    {
        /*
         * Book is "Numerical Analysis 6th Edition by R.L. Burden and J. D. Faires
         * This example is on P291 Table 5.9
         */
        RkfTestOde ode;

        double max_step_size = 0.25;
        double start_time = 0.0;
        double end_time = 2.0;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;

        OdeSolution solutions;

        std::vector<double> state_variables = ode.GetInitialConditions();
        double tolerance = 1e-5;
        solutions = rkf_solver.Solve(&ode, state_variables, start_time, end_time, max_step_size, tolerance);

        // Times (from MatLab Code) to check timstepping is being adapted properly
        TS_ASSERT_DELTA(solutions.rGetTimes()[0], 0, 1e-7);
        TS_ASSERT_DELTA(solutions.rGetTimes()[1], 2.500000000000000e-01, 1e-7);
        TS_ASSERT_DELTA(solutions.rGetTimes()[2], 4.868046415733731e-01, 1e-7);
        TS_ASSERT_DELTA(solutions.rGetTimes()[3], 7.298511818781566e-01, 1e-7);
        TS_ASSERT_DELTA(solutions.rGetTimes()[4], 9.798511818781566e-01, 1e-7);
        TS_ASSERT_DELTA(solutions.rGetTimes()[5], 1.229851181878157e+00, 1e-7);
        TS_ASSERT_DELTA(solutions.rGetTimes()[6], 1.479851181878157e+00, 1e-7);
        TS_ASSERT_DELTA(solutions.rGetTimes()[7], 1.729851181878157e+00, 1e-7);
        TS_ASSERT_DELTA(solutions.rGetTimes()[8], 1.979851181878157e+00, 1e-7);
        TS_ASSERT_DELTA(solutions.rGetTimes()[9], 2.000000000000000e+00, 1e-7);

        TS_ASSERT_EQUALS(solutions.GetNumberOfTimeSteps(), 9u);

        // y values (from analytic result)
        for (unsigned i=0; i<solutions.GetNumberOfTimeSteps(); i++)
        {
            double time = solutions.rGetTimes()[i];
            double y = (time+1.0)*(time+1.0) - 0.5*exp(time);
            // Tolerance set to 1e-5, so 2e-5 to pass here...
            TS_ASSERT_DELTA(solutions.rGetSolutions()[i][0], y, 2e-5);
        }
    }

    void TestRKFehlbergSystemOf3Equations() throw(Exception)
    {
        OdeThirdOrder ode_system;

        double h_value=0.1;

        //Euler solver solution worked out
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;

        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = rkf_solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.25, 1e-5);
        unsigned last = solutions.GetNumberOfTimeSteps();
        double numerical_solution[3];
        numerical_solution[0] = solutions.rGetSolutions()[last][0];
        numerical_solution[1] = solutions.rGetSolutions()[last][1];
        numerical_solution[2] = solutions.rGetSolutions()[last][2];

        // The tests
        double analytical_solution[3];
        analytical_solution[0] = -sin(2);
        analytical_solution[1] = sin(2)+cos(2);
        analytical_solution[2] = 2*sin(2);
        double global_error_rkf = 0.5*2*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(numerical_solution[0],analytical_solution[0],global_error_rkf);
        TS_ASSERT_DELTA(numerical_solution[1],analytical_solution[1],global_error_rkf);
        TS_ASSERT_DELTA(numerical_solution[2],analytical_solution[2],global_error_rkf);
    }

    void TestRKFehlbergNonlinearEquation() throw(Exception)
    {
        Ode4 ode_system;

        double h_value=0.1;

        //Euler solver solution worked out
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();

        solutions = rkf_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, 1e-5);
        int last = solutions.GetNumberOfTimeSteps();

        double numerical_solution;
        numerical_solution = solutions.rGetSolutions()[last][0];

        // The tests
        double analytical_solution = 1.0/(1.0+exp(-12.5));

        TS_ASSERT_DELTA(numerical_solution,analytical_solution,1.0e-4);
    }

    void TestRKFehlbergAnotherNonlinearEquation() throw(Exception)
    {
        Ode5 ode_system;

        double h_value=0.1;
        double end_time = 1.0;

        //Euler solver solution worked out
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();

        solutions = rkf_solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, 1e-5);
        int last = solutions.GetNumberOfTimeSteps();

        double numerical_solution;
        numerical_solution = solutions.rGetSolutions()[last][0];

        // The tests
        double analytical_solution = 1.0/(1.0+4.0*exp(-100.0*end_time));

        TS_ASSERT_DELTA(numerical_solution,analytical_solution,1.0e-3);
    }

    void TestRKFehlbergSystemOf3EquationsWithEvents() throw(Exception)
    {
        OdeThirdOrderWithEvents ode_system_with_events;

        double h_value=0.1;

        //Euler solver solution worked out
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system_with_events.GetInitialConditions();
        solutions = rkf_solver.Solve(&ode_system_with_events, state_variables, 0.0, 2.0, h_value, 1e-5);
        unsigned last = solutions.GetNumberOfTimeSteps();

//        for (unsigned i=0; i<last+1; i++)
//        {
//            std::cout << "Time = " << solutions.rGetTimes()[i] <<
//                " x = " << solutions.rGetSolutions()[i][0] << "\n" << std::flush;
//        }

        // final time should be pi/6 (?)
        TS_ASSERT_DELTA( solutions.rGetTimes()[last], M_PI/6.0, h_value);

        // penultimate y0 should be greater than -0.5
        TS_ASSERT_LESS_THAN(-0.5,solutions.rGetSolutions()[last-1][0]);
        // final y0 should be less than -0.5
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[last][0], -0.5);

        // solver should correctly state the stopping event occured
        TS_ASSERT_EQUALS(rkf_solver.StoppingEventOccured(), true);

        // Coverage of exceptions
        TS_ASSERT_THROWS_ANYTHING(rkf_solver.Solve(&ode_system_with_events, solutions.rGetSolutions()[last], M_PI/6.0, 2.0, 0.1, 1e-5));
        TS_ASSERT_THROWS_ANYTHING(rkf_solver.Solve(&ode_system_with_events, solutions.rGetSolutions()[last],M_PI/6.0, 2.0, 0.1));
    }

    void TestRKFehlbergAnotherNonlinearEquationAnalytic() throw(Exception)
    {
        Ode5Jacobian ode_system;

        double h_value = 0.1;
        double end_time = 1.0;

        //Euler solver solution worked out
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();

        solutions = rkf_solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, 1e-5);
        unsigned last = solutions.GetNumberOfTimeSteps();

        double numerical_solution;
        numerical_solution = solutions.rGetSolutions()[last][0];

        // The tests
        double analytical_solution = 1.0/(1.0+4.0*exp(-100.0*end_time));

        TS_ASSERT_DELTA(numerical_solution,analytical_solution,1.0e-3);
    }

    void TestRKFehlbergVanDerPolOde() throw(Exception)
    {
        VanDerPolOde ode_system;

        double h_value = 1.0;
        double end_time = 100.0;

        //Euler solver solution worked out
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        ode_system.rGetStateVariables() = state_variables;

        rkf_solver.SolveAndUpdateStateVariable(&ode_system, 0.0, end_time, h_value);

        std::vector<double> numerical_solution;
        numerical_solution.reserve(ode_system.GetNumberOfStateVariables());
        numerical_solution = ode_system.rGetStateVariables();

//        OutputFileHandler handler("");
//        out_stream rabbit_file=handler.OpenOutputFile("foxrabbit.dat");
//
//        for (unsigned i=0; i<last; i++)
//        {
//            (*rabbit_file) << solutions.rGetSolutions()[i][0] << "\t" << solutions.rGetSolutions()[i][1] << "\n" << std::flush;
//        }
//        rabbit_file->close();

        // assert that we are within a [-2,2] in x and [-2,2] in y (on limit cycle)
        TS_ASSERT_DELTA(numerical_solution[0],0,2);
        TS_ASSERT_DELTA(numerical_solution[1],0,2);


    }


};

#endif /*TESTRUNGEKUTTAFEHLBERGIVPODESOLVER_HPP_*/
