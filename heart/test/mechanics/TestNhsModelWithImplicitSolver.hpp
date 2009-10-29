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


#ifndef TESTNHSMODELWITHIMPLICITSOLVER_HPP_
#define TESTNHSMODELWITHIMPLICITSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "NhsModelWithImplicitSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ZeroStimulus.hpp"


class TestNhsModelWithImplicitSolver : public CxxTest::TestSuite
{
private:
    double GetSampleCaIValue()
    {
        boost::shared_ptr<EulerIvpOdeSolver> p_euler_solver(new EulerIvpOdeSolver);
        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus);
        LuoRudyIModel1991OdeSystem lr91(p_euler_solver, p_zero_stimulus);
        return lr91.rGetStateVariables()[lr91.GetStateVariableNumberByName("CaI")];
    }

public:
    void TestSolverSingleTimestep()
    {
        NhsModelWithImplicitSolver system_with_solver;

        // lam=const, dlamdt not zero doesn't make much sense, just for testing purposes
        system_with_solver.SetStretchAndStretchRate(0.5, 0.1);
        double Ca_I = GetSampleCaIValue();
        system_with_solver.SetIntracellularCalciumConcentration(Ca_I);

        // solve system (but don't update state vars yet
        system_with_solver.RunDoNotUpdate(0, 0.1, 0.1); // one timestep

        NhsContractionModel system_for_euler_solver;

        unsigned num_vars = system_with_solver.GetNumberOfStateVariables();
        for(unsigned i=0; i<num_vars; i++)
        {
            // both should be the same (ie initial values)
            TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],
                            system_for_euler_solver.rGetStateVariables()[i],
                            1e-12);
        }

        // solve system with euler
        system_for_euler_solver.SetStretchAndStretchRate(0.5, 0.1);
        system_for_euler_solver.SetIntracellularCalciumConcentration(Ca_I);
        EulerIvpOdeSolver euler_solver;
        euler_solver.SolveAndUpdateStateVariable(&system_for_euler_solver, 0, 0.1, 0.1);  // one timestep

        // update state vars on implicit system
        system_with_solver.UpdateStateVariables();

        for(unsigned i=0; i<num_vars; i++)
        {
            //std::cout << system_with_solver.rGetStateVariables()[i] << " "
            //          << system_for_euler_solver.rGetStateVariables()[i] << "\n";

            // we want these within 10% of each other. Note we expect the implicit
            // solver to be more accurate than the explicit solver, and the timestep
            // is quite large (as we want non-zero solutions), so can't expect them
            // to be too close.
            TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],
                            system_for_euler_solver.rGetStateVariables()[i],
                            fabs(system_with_solver.rGetStateVariables()[i]*1e-1));
        }
    }


    void TestSolverManyTimestepsCompareWithEuler()
    {
        for(unsigned run=0; run<2; run++)
        {
            clock_t ck_start, ck_end;

            NhsModelWithImplicitSolver system_with_solver;

            // lam=const, dlamdt not zero doesn't make much sense, just for testing purposes
            system_with_solver.SetStretchAndStretchRate(0.5, 0.1);
            double Ca_I = GetSampleCaIValue();
            // bigger Ca_I, so we get some active tension (and so the a few iterations are
            // needed when solving for T_a and z
            system_with_solver.SetIntracellularCalciumConcentration(10*Ca_I);

            // IF THE SECOND RUN, use the implicit-explicit mixed method for z
            if(run==1)
            {
                system_with_solver.UseImplicitExplicitSolveForZ();
            }

            // solve system and update
            ck_start = clock();
            system_with_solver.RunDoNotUpdate(0, 100, 0.01);
            ck_end = clock();
            double implicit_solve_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

            system_with_solver.UpdateStateVariables();

            // GetNextActiveTension should now be equal to baseclass::GetActiveTension(),
            // as the state vars have been updated
            TS_ASSERT_DELTA(system_with_solver.GetNextActiveTension(),
                            system_with_solver.GetActiveTension(),
                            1e-12);

            // solve system with euler
            NhsContractionModel system_for_euler_solver;
            system_for_euler_solver.SetStretchAndStretchRate(0.5, 0.1);
            system_for_euler_solver.SetIntracellularCalciumConcentration(10*Ca_I);
            EulerIvpOdeSolver euler_solver;

            ck_start = clock();
            euler_solver.SolveAndUpdateStateVariable(&system_for_euler_solver, 0, 100, 0.01);
            ck_end = clock();
            double explicit_solve_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

            run==0 ? std::cout<<"\nImplicit vs Explicit\n" : std::cout<<"\nImplicitExplicit vs Explicit\n";
            unsigned num_vars = system_with_solver.GetNumberOfStateVariables();
            for(unsigned i=0; i<num_vars; i++)
            {
                std::cout << system_with_solver.rGetStateVariables()[i] << " "
                          << system_for_euler_solver.rGetStateVariables()[i] << "\n";

                // small timesteps, want these to be very close (1%). and they are
                TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],
                                system_for_euler_solver.rGetStateVariables()[i],
                                fabs(system_with_solver.rGetStateVariables()[i]*1e-2));
            }
            std::cout << "TIMES: " << implicit_solve_time << " " << explicit_solve_time << "\n\n";

            // for coverage
            system_with_solver.SetActiveTensionInitialGuess(system_with_solver.GetActiveTension());
        }
    }

//// test how large a timestep the implicit solver can get away with. needs more study
//    void TestImplicitSolverWithLargeTimeSteps()
//    {
//        NhsModelWithImplicitSolver system_with_solver;
//        system_with_solver.SetStretchAndStretchRate(0.5, 0.1);
//        system_with_solver.SetIntracellularCalciumConcentration(10*GetSampleCaIValue());
//
//        system_with_solver.RunDoNotUpdate(0, 100, 0.01);
//        system_with_solver.UpdateStateVariables();
//
//        NhsModelWithImplicitSolver system_with_solver2;
//        system_with_solver2.SetStretchAndStretchRate(0.5, 0.1);
//        system_with_solver2.SetIntracellularCalciumConcentration(10*GetSampleCaIValue());
//
//        system_with_solver2.RunDoNotUpdate(0, 100, 1);
//        system_with_solver2.UpdateStateVariables();
//
//        unsigned num_vars = system_with_solver.GetNumberOfStateVariables();
//        for(unsigned i=0; i<num_vars; i++)
//        {
//            std::cout << system_with_solver.rGetStateVariables()[i] << " "
//                      << system_with_solver2.rGetStateVariables()[i] << "\n";
//
//            TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],
//                            system_with_solver2.rGetStateVariables()[i],
//                            fabs(system_with_solver.rGetStateVariables()[i]*1e-2));
//        }
//    }

    // test that checks RunDoNotUpdate does not do anything permanent on the class,
    // by checking by doing things repeatedly, and changing the order, makes no
    // difference
    void TestRunDoesNotUpdate()
    {
        NhsModelWithImplicitSolver system;

        double Ca_I = GetSampleCaIValue();
        system.SetIntracellularCalciumConcentration(Ca_I);

        // get initial active tension
        double init_Ta = system.GetActiveTension();

        system.SetStretchAndStretchRate(0.6, 0.1);
        system.RunDoNotUpdate(0, 1, 0.01);

        double Ta1 = system.GetNextActiveTension();

        system.SetStretchAndStretchRate(0.6, 0.2);
        system.RunDoNotUpdate(0, 1, 0.01);

        double Ta2 = system.GetNextActiveTension();

        // note that lam/end time etc must be large enough for there
        // to be non-zero Ta at the next time
        TS_ASSERT_DIFFERS(init_Ta, Ta1);
        TS_ASSERT_DIFFERS(init_Ta, Ta2);

        system.SetStretchAndStretchRate(0.6, 0.2);
        system.RunDoNotUpdate(0, 1, 0.01);

        double should_be_Ta2 = system.GetNextActiveTension();

        system.SetStretchAndStretchRate(0.6, 0.1);
        system.RunDoNotUpdate(0, 1, 0.01);

        double should_be_Ta1 = system.GetNextActiveTension();

        TS_ASSERT_EQUALS(Ta1, should_be_Ta1);
        TS_ASSERT_EQUALS(Ta2, should_be_Ta2);

        system.SetStretchAndStretchRate(0.6, 0.1);
        system.RunDoNotUpdate(0, 1, 0.01);

        double should_also_be_Ta1 = system.GetNextActiveTension();

        system.SetStretchAndStretchRate(0.6, 0.2);
        system.RunDoNotUpdate(0, 1, 0.01);

        double should_also_be_Ta2 = system.GetNextActiveTension();

        TS_ASSERT_EQUALS(Ta1, should_also_be_Ta1);
        TS_ASSERT_EQUALS(Ta2, should_also_be_Ta2);
    }

    void TestGetActiveTension()
    {
        NhsModelWithImplicitSolver system;

        double Ca_I = GetSampleCaIValue();
        system.SetIntracellularCalciumConcentration(Ca_I);
        system.SetStretchAndStretchRate(0.6, 0.1);
        system.RunDoNotUpdate(0, 1, 0.01);

        double Ta_at_next_time_before_update = system.GetNextActiveTension();

        system.UpdateStateVariables();

        double Ta = system.GetActiveTension();

        TS_ASSERT_DELTA(Ta, Ta_at_next_time_before_update, 1e-12);
    }
};

#endif /*TESTNHSMODELWITHIMPLICITSOLVER_HPP_*/
