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
#ifndef TESTTYSONNOVAK2001ODESYSTEM_HPP_
#define TESTTYSONNOVAK2001ODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <stdio.h>
#include <ctime>
#include <vector>
#include <iostream>

#include "TysonNovak2001OdeSystem.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestTysonNovak2001OdeSystem : public CxxTest::TestSuite
{
public:

    void TestTysonNovakEquation()
    {
        TysonNovak2001OdeSystem tyson_novak_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(0.6);
        initial_conditions.push_back(0.1);
        initial_conditions.push_back(1.5);
        initial_conditions.push_back(0.6);
        initial_conditions.push_back(0.6);
        initial_conditions.push_back(0.85);

        std::vector<double> derivs(initial_conditions.size());
        tyson_novak_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct
        // Divided by 60 to change to hours
        TS_ASSERT_DELTA(derivs[0],-4.400000000000000e-02*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-6.047872340425530e+00*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2],3.361442884485838e-02*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[3],4.016602000735009e-02*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[4],8.400000000000001e-03*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[5],7.777500000000001e-03*60.0, 1e-5);
    }

    void TestTysonNovakSolver() throw(Exception)
    {
        TysonNovak2001OdeSystem tyson_novak_system;

        // Solve system using backward Euler solver

        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits

        double dt=0.1/60.0;

        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(6);

        OdeSolution solutions;

        std::vector<double> state_variables = tyson_novak_system.GetInitialConditions();

        double start_time, end_time, elapsed_time = 0.0;

        state_variables = tyson_novak_system.GetInitialConditions();

        start_time = std::clock();
        solutions = backward_euler_solver.Solve(&tyson_novak_system, state_variables, 0.0, 75.8350/60.0, dt, dt);
        end_time = std::clock();

        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Elapsed time = " << elapsed_time << "\n";

        // If you run it up to about 75min the ODE will stop, anything less and it will not and this test will fail
        TS_ASSERT(backward_euler_solver.StoppingEventOccurred());

        unsigned end = solutions.rGetSolutions().size() - 1;

        // The following code provides nice output for gnuplot
        // use the command
        // plot "tyson_novak.dat" u 1:2
        // or
        // plot "tyson_novak.dat" u 1:3 etc. for the various proteins...

//        OutputFileHandler handler("");
//        out_stream file=handler.OpenOutputFile("tyson_novak.dat");
//        for (unsigned i=0; i<=end; i++)
//        {
//            (*file) << solutions.rGetTimes()[i]<< "\t" << solutions.rGetSolutions()[i][0] << "\t" << solutions.rGetSolutions()[i][1] <<"\t"<< solutions.rGetSolutions()[i][2] << "\t" << solutions.rGetSolutions()[i][3] <<"\t"<< solutions.rGetSolutions()[i][4] << "\t" << solutions.rGetSolutions()[i][5] << "\n" << std::flush;
//        }
//        file->close();

        int my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if (my_rank==0) // if master process
        {

            int step_per_row = 1;
            ColumnDataWriter writer("TysonNovak","TysonNovak");
            int time_var_id = writer.DefineUnlimitedDimension("Time","s");

            std::vector<int> var_ids;
            for (unsigned i=0; i<tyson_novak_system.rGetVariableNames().size(); i++)
            {
                var_ids.push_back(writer.DefineVariable(tyson_novak_system.rGetVariableNames()[i],
                                                        tyson_novak_system.rGetVariableUnits()[i]));
            }
            writer.EndDefineMode();

            for (unsigned i = 0; i < solutions.rGetSolutions().size(); i+=step_per_row)
            {
                writer.PutVariable(time_var_id, solutions.rGetTimes()[i]);
                for (unsigned j=0; j<var_ids.size(); j++)
                {
                    writer.PutVariable(var_ids[j], solutions.rGetSolutions()[i][j]);
                }
                writer.AdvanceAlongUnlimitedDimension();
            }
            writer.Close();
        }
        MPI_Barrier(PETSC_COMM_WORLD);

        // Proper values calculated using the Matlab stiff ODE solver ode15s. Note that
        // large tolerances are required for the tests to pass with both chaste solvers
        // and CVODE.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],1.54216806705641, 1e-1);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.40562614481544, 1e-1);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],0.95328206604519, 2e-2);

    }
};

#endif /*TESTTYSONNOVAK2001ODESYSTEM_HPP_*/
