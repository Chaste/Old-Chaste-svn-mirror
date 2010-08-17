/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef TESTMIRAMS2010WNTODESYSTEM_HPP_
#define TESTMIRAMS2010WNTODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>
#include <fstream>
#include <ctime>

#include "OutputFileHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "Mirams2010WntOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

class TestMirams2010WntOdeSystem : public AbstractCellBasedTestSuite
{
public:

    void TestMirams2010WntOdeSystemSetup() throw(Exception)
    {
#ifdef CHASTE_CVODE
        double wnt_level = 0.5;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        Mirams2010WntOdeSystem wnt_system(wnt_level, p_state);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value = 0.0001;

        CvodeAdaptor cvode_solver;

        OdeSolution solutions;
        //OdeSolution solutions2;

        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        double start_time, end_time, elapsed_time = 0.0;
        start_time = std::clock();
        solutions = cvode_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Cvode Elapsed time = " << elapsed_time << " secs for 100 hours\n";

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100, 1e-2);

        // Decent results
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 67.5011, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 67.5011, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], wnt_level, 1e-4);
#endif //CHASTE_CVODE
    }

    void TestGarysWntOdeSystemApc2Hit() throw(Exception)
    {
#ifdef CHASTE_CVODE
        double wnt_level = 0.5;
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        Mirams2010WntOdeSystem wnt_system(wnt_level, p_apc2);

        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits

        double h_value = 0.0001;
        CvodeAdaptor cvode_solver;
        OdeSolution solutions;
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        double start_time, end_time, elapsed_time = 0.0;
        start_time = std::clock();
        solutions = cvode_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Cvode Elapsed time = " << elapsed_time << " secs for 100 hours\n";

        // Test solutions are OK for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Test the simulation is ending at the right time (going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100, 1e-2);

        // Check results are correct
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 433.1155, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 433.1155, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], wnt_level, 1e-4);
#endif //CHASTE_CVODE
    }

    void TestGarysWntOdeSystemBetaCatenin1Hit() throw(Exception)
    {
#ifdef CHASTE_CVODE
        double wnt_level = 0.5;
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        Mirams2010WntOdeSystem wnt_system(wnt_level, p_bcat1);

        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits

        double h_value = 0.0001;
        CvodeAdaptor cvode_solver;
        OdeSolution solutions;
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();

        double start_time, end_time, elapsed_time = 0.0;
        start_time = std::clock();
        solutions = cvode_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Cvode Elapsed time = " << elapsed_time << " secs for 100 hours\n";

        // Test solutions are OK for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Tests the simulation is ending at the right time (going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100, 1e-2);

        // Check results are correct
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 67.5011, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 824.0259, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], wnt_level, 1e-4);
#endif //CHASTE_CVODE
    }
};

#endif /* TESTMIRAMS2010WNTODESYSTEM_HPP_ */
