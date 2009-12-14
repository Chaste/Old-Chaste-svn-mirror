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


#ifndef _TESTIONICMODELS_HPP_
#define _TESTIONICMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "RunAndCheckIonicModels.hpp"
#include "Exception.hpp"

#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"
#include "MultiStimulus.hpp"
#include "ZeroStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CellProperties.hpp"

#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BackwardEulerLuoRudyIModel1991.hpp"

#include "FoxModel2002Modified.hpp"
#include "BackwardEulerFoxModel2002Modified.hpp"

#include "FaberRudy2000Version3.hpp"
#include "FaberRudy2000Version3Optimised.hpp"

#include "NobleVargheseKohlNoble1998.hpp"
#include "NobleVargheseKohlNoble1998Optimised.hpp"
#include "BackwardEulerNobleVargheseKohlNoble1998.hpp"
#include "Mahajan2008OdeSystem.hpp"
#include "TenTusscher2006OdeSystem.hpp"
#include "DiFrancescoNoble1985OdeSystem.hpp"
#include "Maleckar2009OdeSystem.hpp"
#include "ArchiveLocationInfo.hpp"

#include "PetscTools.hpp" //No PETSc here -- this is just to double-check

// Note: RunOdeSolverWithIonicModel(), CheckCellModelResults(), CompareCellModelResults()
// are defined in RunAndCheckIonicModels.hpp

class TestIonicModels : public CxxTest::TestSuite
{
public:
    void TestOdeSolveForNoble98WithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -3;  // uA/cm2
        double duration_stimulus = 3;  // ms
        double start_stimulus = 10.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(
                magnitude_stimulus,
                duration_stimulus,
                start_stimulus));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        double time_step = 0.01;

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

        //Check Standard
        CML_noble_varghese_kohl_noble_1998_basic n98_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_ode_system,
                                   150.0,
                                   "N98RegResult");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("N98RegResult");
        TS_ASSERT_DELTA( n98_ode_system.GetIIonic(), 0.2462, 1e-3);

        std::string error_should_be = "Non fast-slow cell model being used in a fast-slow problem.";
        TS_ASSERT_THROWS_THIS(n98_ode_system.SetState(ALL_VARS), error_should_be);
        TS_ASSERT_THROWS_THIS(n98_ode_system.GetNumSlowValues(), error_should_be);
        TS_ASSERT_THROWS_THIS(n98_ode_system.IsFastOnly(), error_should_be);
        std::vector<double> slows;
        slows.push_back(100.0);
        TS_ASSERT_THROWS_THIS(n98_ode_system.AdjustOutOfRangeSlowValues(slows), error_should_be);
        TS_ASSERT_THROWS_THIS(n98_ode_system.GetSlowValues(slows), error_should_be);
        TS_ASSERT_THROWS_THIS(n98_ode_system.SetSlowValues(slows), error_should_be);

     }

    void TestOdeSolveForOptimisedNoble98WithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -3;  // uA/cm2
        double duration_stimulus = 3;  // ms
        double start_stimulus = 10.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                        duration_stimulus,
                                                                        start_stimulus));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        double time_step = 0.01;

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

        //Check Optimised
        CML_noble_varghese_kohl_noble_1998_basic_pe_lut n98_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_ode_system,
                                   150.0,
                                   "N98RegResultOpt");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("N98RegResultOpt", "N98RegResult");
        TS_ASSERT_DELTA( n98_ode_system.GetIIonic(), 0.2462, 1e-3);

        //Stress the lookup table with a silly voltage
        n98_ode_system.rGetStateVariables()[0] = 70.0;
        TS_ASSERT_EQUALS(n98_ode_system.GetVoltage(), 70.0);
        TS_ASSERT_THROWS_EQUALS( n98_ode_system.GetIIonic(), const Exception &err,
                err.GetShortMessage().find("V outside lookup table range",0), 0u);
        n98_ode_system.rGetStateVariables()[0] = 71.0;
        TS_ASSERT_THROWS_EQUALS( n98_ode_system.GetIIonic(), const Exception &err,
                err.GetShortMessage().find("V outside lookup table range",0), 0u);
        n98_ode_system.rGetStateVariables()[0] = 69.0;
        TS_ASSERT_THROWS_NOTHING( n98_ode_system.GetIIonic());
        n98_ode_system.rGetStateVariables()[0] = -100.1;
        TS_ASSERT_THROWS_EQUALS( n98_ode_system.GetIIonic(), const Exception &err,
                err.GetShortMessage().find("V outside lookup table range",0), 0u);
        n98_ode_system.rGetStateVariables()[0] = -100.0;
        TS_ASSERT_THROWS_NOTHING( n98_ode_system.GetIIonic());

        n98_ode_system.rGetStateVariables()[0] = -100.1;

        TS_ASSERT_THROWS_EQUALS( RunOdeSolverWithIonicModel(&n98_ode_system, 150.0, "DoNotRun"),
                const Exception &err, err.GetShortMessage().find("V outside lookup table range",0), 0u);
    }


    void TestOdeSolverForHH52WithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = 20.0;  // uA/cm2
        double duration_stimulus = 0.5;  // ms
        double start_stimulus = 10.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                        duration_stimulus,
                                                                        start_stimulus));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&hh52_ode_system,
                                   150.0,
                                   "HH52RegResult");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("HH52RegResult");

        // test GetIionic: (the GetIionic method was first manually tested
        // by changing the EvaluateYDerivatives() code to call it, this verified
        // that GetIionic has no errors, therefore we can test here against
        // a hardcoded result
        RunOdeSolverWithIonicModel(&hh52_ode_system,
                                   15.0,
                                   "HhGetIIonic");
        TS_ASSERT_DELTA( hh52_ode_system.GetIIonic(), 40.6341, 1e-3);
    }


    void TestOdeSolverForFHN61WithSimpleStimulus(void) throw (Exception)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -80.0;   // dimensionless
        double duration_stimulus = 0.5;  // ms
        double start_stimulus = 0.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                        duration_stimulus,
                                                                        start_stimulus));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        FitzHughNagumo1961OdeSystem fhn61_ode_system(p_solver, p_stimulus);

        // fhn has no [Ca_i]
        TS_ASSERT_THROWS_THIS(fhn61_ode_system.GetIntracellularCalciumConcentration(),
                "AbstractCardiacCell::GetIntracellularCalciumConcentration() called. "
                "Either model has no [Ca_i] or method has not been implemented yet");

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&fhn61_ode_system,
                                   500.0,
                                   "FHN61RegResult");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("FHN61RegResult");

        // test GetIionic ('fake' ionic current) (the GetIionic method was first
        // manually tested by changing the EvaluateYDerivatives() code to call it,
        // this verified that GetIionic has no errors, therefore we can test here
        // against a hardcoded result
        TS_ASSERT_DELTA( fhn61_ode_system.GetIIonic(), -0.0058, 1e-3);

        // some coverage
        boost::shared_ptr<SimpleStimulus> p_another_stimulus(new SimpleStimulus(-200,1.0, 0.0));
        boost::shared_ptr<SimpleStimulus> p_intra_stimulus(new SimpleStimulus(-100,1.0, 0.0));
        FitzHughNagumo1961OdeSystem another_fhn61_ode_system(p_solver, p_stimulus);

        another_fhn61_ode_system.SetStimulusFunction(p_another_stimulus);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetStimulus(0.5), -200, 1e-12);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetIntracellularStimulus(0.5), -200, 1e-12);

        another_fhn61_ode_system.SetIntracellularStimulusFunction(p_intra_stimulus);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetStimulus(0.5), -100, 1e-12);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetIntracellularStimulus(0.5), -100, 1e-12);
    }


    void TestOdeSolverForLR91WithDelayedSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        double end_time = 1000.0; //One second in milliseconds

        LuoRudyIModel1991OdeSystem lr91_ode_system(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(lr91_ode_system.GetVoltageIndex(), 4u); // For coverage

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91DelayedStim");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;

        CheckCellModelResults("Lr91DelayedStim");

        // test GetIionic: (the GetIionic method was first manually tested
        // by changing the EvaluateYDerivatives() code to call it, this verified
        // that GetIionic has no errors, therefore we can test here against
        // a hardcoded result
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   60.0,
                                   "Lr91GetIIonic");
        TS_ASSERT_DELTA( lr91_ode_system.GetIIonic(), 1.9411, 1e-3);
    }

    void TestOdeSolverForLR91WithRegularStimulus(void) throw (Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude,
                                                                          duration,
                                                                          period,
                                                                          start));

        double end_time = 1000.0; //One second in milliseconds

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        LuoRudyIModel1991OdeSystem lr91_ode_system(p_solver, p_stimulus);

        // cover get intracellular calcium
        TS_ASSERT_DELTA(lr91_ode_system.GetIntracellularCalciumConcentration(), 0.0002, 1e-5)

        // Solve and write to file
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91RegularStim");

        CheckCellModelResults("Lr91RegularStim");
    }

    void TestBackwardEulerLr91WithDelayedSimpleStimulus(void) throw (Exception)
    {
        clock_t ck_start, ck_end;
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        double end_time = 1000.0; //One second in milliseconds

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Solve using backward euler
        BackwardEulerLuoRudyIModel1991 lr91_backward_euler(p_stimulus);

        // cover get intracellular calcium
        TS_ASSERT_DELTA(lr91_backward_euler.GetIntracellularCalciumConcentration(), 0.0002, 1e-5)

        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_backward_euler,
                                   end_time,
                                   "Lr91BackwardEuler");
        ck_end = clock();
        double backward1 = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        // Solve using forward Euler
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        LuoRudyIModel1991OdeSystem lr91_ode_system(p_solver, p_stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91DelayedStim");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        // Compare results
        CompareCellModelResults("Lr91DelayedStim", "Lr91BackwardEuler", 0.01);

        // Try with larger timestep and coarser tolerance.
        // We can't use a larger time step than 0.01 for forward Euler - the gating
        // variables go out of range.

        // (Use alternative contructor for coverage. This is a hack -see ticket:451 )

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.5, 0.5, 0.5);
        BackwardEulerLuoRudyIModel1991 lr91_backward_euler2(p_solver, p_stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_backward_euler2,
                                   end_time,
                                   "Lr91BackwardEuler2");

        ck_end = clock();
        double backward2 = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        CompareCellModelResults("Lr91DelayedStim", "Lr91BackwardEuler2", 0.25);

        std::cout << "Run times:\n\tForward: " << forward << "\n\tBackward: "
                  << backward1 << "\n\tBackward (long dt): " << backward2 << std::endl;


        // cover and check GetIIonic() match for normal and backward euler lr91
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);
        LuoRudyIModel1991OdeSystem lr91(p_solver, p_stimulus);
        BackwardEulerLuoRudyIModel1991 backward_lr91(p_stimulus);
        // calc IIonic using initial conditions
        TS_ASSERT_DELTA(lr91.GetIIonic(), backward_lr91.GetIIonic(), 1e-12);

        // cover alternative constructor
        BackwardEulerLuoRudyIModel1991 lr91_backward_euler3(p_solver, p_stimulus);
    }

    void TestOdeSolverForFR2000WithDelayedSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0;  // ms
        double when = 10.0; // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        double end_time = 1000.0; //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.007, 0.007, 0.007);

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        FaberRudy2000Version3Optimised fr2000_ode_system_opt(p_solver, p_stimulus);
        FaberRudy2000Version3 fr2000_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&fr2000_ode_system,
                                   end_time,
                                   "FR2000DelayedStim",
                                   500, false);
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        ck_start = clock();
        RunOdeSolverWithIonicModel(&fr2000_ode_system_opt,
                                   end_time,
                                   "FR2000DelayedStimOpt",
                                   500, false);
        ck_end = clock();
        double opt = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        std::cout << "\n\tForward: " << forward
                  << "\n\tOptimised: " << opt << std::endl;

        CheckCellModelResults("FR2000DelayedStim");
        CompareCellModelResults("FR2000DelayedStim", "FR2000DelayedStimOpt", 1e-4);

        // test GetIionic: (the GetIionic method was first manually tested
        // by changing the EvaluateYDerivatives() code to call it, this verified
        // that GetIionic has no errors, therefore we can test here against
        // a hardcoded result
        TS_ASSERT_DELTA(fr2000_ode_system.GetIIonic(), 0.0002, 1e-4);
        TS_ASSERT_DELTA(fr2000_ode_system_opt.GetIIonic(), 0.0002, 1e-4);

        //Check that ComputeExceptVoltage does the correct thing (doesn't change the voltage)
        double voltage=fr2000_ode_system.GetVoltage();
        fr2000_ode_system.ComputeExceptVoltage(end_time, end_time+0.001);
        TS_ASSERT_DELTA(fr2000_ode_system.GetVoltage(), voltage, 1e-5);
        voltage=fr2000_ode_system.GetVoltage();
        fr2000_ode_system_opt.ComputeExceptVoltage(end_time, end_time+0.001);
        TS_ASSERT_DELTA(fr2000_ode_system_opt.GetVoltage(), voltage, 1e-5);

    }


    void TestOdeSolverForFR2000WithVariablePotassiumCurrents(void)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0;  // ms
        double when = 0.0; // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        double end_time = 1000.0; //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.007, 0.007, 0.007);

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        FaberRudy2000Version3 fr2000_ode_system_endo(p_solver, p_stimulus);
        fr2000_ode_system_endo.SetScaleFactorGks(0.462);
        fr2000_ode_system_endo.SetScaleFactorIto(0.0);
        fr2000_ode_system_endo.SetScaleFactorGkr(1.0);
        // Solve and write to file
        RunOdeSolverWithIonicModel(&fr2000_ode_system_endo,
                                   end_time,
                                   "FR2000Endo",
                                   500, false);

        CheckCellModelResults("FR2000Endo");

        FaberRudy2000Version3 fr2000_ode_system_mid(p_solver, p_stimulus);
        fr2000_ode_system_mid.SetScaleFactorGks(1.154);
        fr2000_ode_system_mid.SetScaleFactorIto(0.85);
        fr2000_ode_system_mid.SetScaleFactorGkr(1.0);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&fr2000_ode_system_mid,
                                   end_time,
                                   "FR2000Mid",
                                   500, false);

        CheckCellModelResults("FR2000Mid");

        FaberRudy2000Version3 fr2000_ode_system_epi(p_solver, p_stimulus);
        fr2000_ode_system_epi.SetScaleFactorGks(1.154);
        fr2000_ode_system_epi.SetScaleFactorIto(1.0);
        fr2000_ode_system_epi.SetScaleFactorGkr(1.0);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&fr2000_ode_system_epi,
                                   end_time,
                                   "FR2000Epi",
                                   500, false);

        CheckCellModelResults("FR2000Epi");
    }


    void TestOdeSolverForFox2002WithRegularStimulus(void) throw (Exception)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -80.0;
        double duration  = 1.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude,
                                                                          duration,
                                                                          period,
                                                                          start));

        double end_time = 200.0;  // milliseconds

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.002, 0.002, 0.002); // 0.005 leads to NaNs.

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        FoxModel2002Modified fox_ode_system(p_solver, p_stimulus);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);
        BackwardEulerFoxModel2002Modified backward_system(p_stimulus);

        // Mainly for coverage, and to test consistency of GetIIonic
        TS_ASSERT_DELTA(fox_ode_system.GetIIonic(),
                        backward_system.GetIIonic(),
                        1e-6);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&fox_ode_system,
                                   end_time,
                                   "FoxRegularStim",
                                   500);
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        CheckCellModelResults("FoxRegularStim");

        // Solve using Backward Euler
        ck_start = clock();
        RunOdeSolverWithIonicModel(&backward_system,
                                   end_time,
                                   "BackwardFoxRegularStim",
                                   100);
        ck_end = clock();
        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        CompareCellModelResults("FoxRegularStim", "BackwardFoxRegularStim", 0.2);

        std::cout << "Run times:\n\tForward: " << forward
                  << "\n\tBackward: " << backward
                  << std::endl;

    }

    // For some bizarre reason having the exception specification here and/or in the private method
    // causes the IntelProduction build to fail on this test (unless it is the only test defined, i.e.
    // we x-out the other tests in this file).
    void TestLr91WithVoltageDropVariousTimeStepRatios() //throw (Exception)
    {
        TS_ASSERT_THROWS_CONTAINS(TryTestLr91WithVoltageDrop(1),
                "m gate for fast sodium current has gone out of range. Check model parameters, for example spatial stepsize\n");
        TS_ASSERT_THROWS_CONTAINS(TryTestLr91WithVoltageDrop(2),
                "m gate for fast sodium current has gone out of range. Check model parameters, for example spatial stepsize\n");
        TS_ASSERT_THROWS_CONTAINS(TryTestLr91WithVoltageDrop(3),
                "m gate for fast sodium current has gone out of range. Check model parameters, for example spatial stepsize\n");
        TS_ASSERT_THROWS_NOTHING(TryTestLr91WithVoltageDrop(4));
    }

    void TestOdeSolveForBackwardNoble98WithSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -3;  // uA/cm2
        double duration_stimulus = 3;  // ms
        double start_stimulus = 10.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                        duration_stimulus,
                                                                        start_stimulus));

        // Just adding to check that multi-stim works properly with a cell model.
        boost::shared_ptr<MultiStimulus> p_multi_stim(new MultiStimulus);
        p_multi_stim->AddStimulus(p_stimulus);

        double time_step = 0.2;

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

        BackwardEulerNobleVargheseKohlNoble1998 n98_backward_system(p_multi_stim);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_backward_system,
                                   150.0,
                                   "N98BackwardResult",
                                   1);
        ck_end = clock();
        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tBackward: " << backward << std::endl;

        CheckCellModelResults("N98BackwardResult");
        //TS_ASSERT_DELTA( n98_backward_system.GetIIonic(), 0.023, 1e-3);
        //This cell now returns a current density
        TS_ASSERT_DELTA( n98_backward_system.GetIIonic(), 0.2462, 1e-3);

        ///\todo compare with the forward results?
    }

    void TestOdeSolveForTT06WithSimpleStimulus(void)
    {

        double simulation_end=350;/*end time, in milliseconds for this model*/

        // Set the stimulus, the following values are appropriate for single cell simulations of this model.
        double magnitude = -38.0;   // pA/pF
        double duration = 1.0;  // ms
        double start = 100;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude,
                                                                        duration,
                                                                        start));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);// with Forward Euler, this must be as small as 0.001.
        TenTusscher2006OdeSystem TT_model(p_solver, p_stimulus);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&TT_model,
                                   simulation_end,
                                   "TenTusscher",
                                   1000,
                                   true);
        //Check against validated data
        //These data are considered valid after (visually) checking against output from  CellML code of the model for an epicardial cell
        // and also numerically compared against pycml automatically generated code.
        CheckCellModelResults("TenTusscher");

        //Test the GetIIonic method against one hardcoded value.
        TS_ASSERT_DELTA( TT_model.GetIIonic(), 0.5279, 1e-3);

        //Test the GetIIonic method against one hardcoded value for initial values of voltage
        //(mainly for coverage of different if conditions in sodium channel gates for different voltages)
        TenTusscher2006OdeSystem TT_model_initial(p_solver, p_stimulus);
        TS_ASSERT_DELTA(TT_model_initial.GetIIonic(), 0.0012 , 1e-3);

        //now test the scale factor methods

        TT_model.SetScaleFactorGks(1.0);
        TT_model.SetScaleFactorIto(1.0);
        //run for only 10 ms
        RunOdeSolverWithIonicModel(&TT_model,
                                   10,
                                   "TenTusscher",
                                   1000,
                                   false);
        double i_ionic = TT_model.GetIIonic();

        //now double the scale factors
        TT_model.SetScaleFactorGks(2.0);
        TT_model.SetScaleFactorIto(2.0);
        //run again for only 10 ms
        RunOdeSolverWithIonicModel(&TT_model,
                                   10,
                                   "TenTusscher",
                                   1000,
                                   false);
        double i_ionic_2 = TT_model.GetIIonic();

        //check that the second case gets a smaller i_ionic
        TS_ASSERT_LESS_THAN(i_ionic , i_ionic_2);

        TT_model.SetScaleFactorGkr(0.0);
        //run again for only 10 ms
        RunOdeSolverWithIonicModel(&TT_model,
                                   10,
                                   "TenTusscher",
                                   1000,
                                   false);
        double i_ionic_3 = TT_model.GetIIonic();

         TS_ASSERT_LESS_THAN(i_ionic , i_ionic_3);
     }

    void TestDifrancescoNoble1985(void) throw (Exception)
    {
        // Set stimulus (no stimulus in this case because this cell is self excitatory)
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.01);
        DiFrancescoNoble1985OdeSystem purkinje_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&purkinje_ode_system,
                                   1800,/*end time, in milliseconds for this model*/
                                   "DiFrancescoNoble",
                                   1000);
        //Check against validated data
        //(the valid data have been checked against CellML code of the model known to be valid).
        CheckCellModelResults("DiFrancescoNoble");

         //Test the GetIIonic method against one hardcoded value.
        TS_ASSERT_DELTA(purkinje_ode_system.GetIIonic(), -0.0141, 1e-3);
     }

    void TestMahajan2008(void) throw (Exception)
    {
        // Set stimulus
        double magnitude_stimulus = -1800;
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude_stimulus,
                                                                          0.05,
                                                                          1000,
                                                                          10.0));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);
        Mahajan2008OdeSystem rabbit_ode_system(p_solver, p_stimulus);

        //Test the GetIIonic method against one hardcoded value.
        TS_ASSERT_DELTA(rabbit_ode_system.GetIIonic(), 0.0027, 1e-3);

        // Solve and write to file
        RunOdeSolverWithIonicModel(&rabbit_ode_system,
                                   400,/*end time, in milliseconds for this model*/
                                   "Mahajan2008",
                                   1000);
        // Check against validated data
        // (the code for the mahajan model was generated from a CellML code known to be valid)
        CheckCellModelResults("Mahajan2008");
    }


    /**
     *  Here we test the scale factors methiods for the mahajan model.
     *  The idea is to set the scale factors for the 3 different cell types (epi, mid and endo)
     *  and check that the rsulting APD makes sense if compared to experiemntal results
     */
    void TestScaleFactorsForMahajanModel(void) throw(Exception)
    {
        double end_time=300;
        double time_step=0.01;
        double sampling_time=time_step;

        // Set stimulus
        double magnitude_stimulus = -70.0;   // pA/pF
        double duration_stimulus = 1.0;  // ms
        double start_stimulus = 10.0;   // ms
        double period=1000;//here, this is ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude_stimulus,
                                                                          duration_stimulus,
                                                                          period,
                                                                          start_stimulus));

        boost::shared_ptr<EulerIvpOdeSolver> p_forward_solver(new EulerIvpOdeSolver); //define the solver
        Mahajan2008OdeSystem forward_model(p_forward_solver, p_stimulus);
        boost::shared_ptr<BackwardEulerIvpOdeSolver> p_backward_solver(new BackwardEulerIvpOdeSolver(
                                    forward_model.GetNumberOfStateVariables()));

        Mahajan2008OdeSystem epicardial_model(p_backward_solver, p_stimulus);
        Mahajan2008OdeSystem endocardial_model(p_backward_solver, p_stimulus);
        Mahajan2008OdeSystem midmyocardial_model(p_backward_solver, p_stimulus);

        epicardial_model.SetScaleFactorGks(1.0);
        epicardial_model.SetScaleFactorIto(1.0);
        epicardial_model.SetScaleFactorGkr(1.0);

        midmyocardial_model.SetScaleFactorGks(0.09);
        midmyocardial_model.SetScaleFactorIto(1.0);
        midmyocardial_model.SetScaleFactorGkr(1.0);

        endocardial_model.SetScaleFactorGks(0.86);
        endocardial_model.SetScaleFactorIto(0.2);
        endocardial_model.SetScaleFactorGkr(1.0);

        std::vector<double> state_variables_epi = epicardial_model.GetInitialConditions();
        std::vector<double> state_variables_endo = endocardial_model.GetInitialConditions();
        std::vector<double> state_variables_mid = midmyocardial_model.GetInitialConditions();

        const std::string mahajan_epi_file = "Mahajan_epi";
        const std::string mahajan_mid_file = "Mahajan_mid";
        const std::string mahajan_endo_file = "Mahajan_endo";

        // Solve and write to file

        OdeSolution epi_solution;
        epi_solution = p_backward_solver->Solve(&epicardial_model, state_variables_epi, 0, end_time, time_step, sampling_time);

        epi_solution.WriteToFile("TestIonicModels",
                                 mahajan_epi_file,
                                 &epicardial_model,
                                 "ms",//time units
                                 100,//steps per row
                                 false);/*true cleans the directory*/

        OdeSolution mid_solution;
        mid_solution = p_backward_solver->Solve(&midmyocardial_model, state_variables_mid, 0, end_time, time_step, sampling_time);

        mid_solution.WriteToFile("TestIonicModels",
                                 mahajan_mid_file,
                                 &epicardial_model,
                                 "ms",//time units
                                 100,//steps per row
                                 false);/*true cleans the directory*/

        OdeSolution endo_solution;
        endo_solution = p_backward_solver->Solve(&endocardial_model, state_variables_endo, 0, end_time, time_step, sampling_time);

        endo_solution.WriteToFile("TestIonicModels",
                                  mahajan_endo_file,
                                  &midmyocardial_model,
                                  "ms",//time units
                                  100,//steps per row
                                  false);/*true cleans the directory*/


        ColumnDataReader data_reader_epi("TestIonicModels", mahajan_epi_file);
        ColumnDataReader data_reader_mid("TestIonicModels", mahajan_mid_file);
        ColumnDataReader data_reader_endo("TestIonicModels", mahajan_endo_file);

        std::vector<double> times = data_reader_epi.GetValues("Time");
        std::vector<double> v_endo = data_reader_endo.GetValues("V");
        std::vector<double> v_epi = data_reader_epi.GetValues("V");
        std::vector<double> v_mid = data_reader_mid.GetValues("V");

        CellProperties  cell_properties_endo(v_endo, times);
        CellProperties  cell_properties_epi(v_epi, times);
        CellProperties  cell_properties_mid(v_mid, times);

        double epi_APD = cell_properties_epi.GetLastActionPotentialDuration(90);
        double endo_APD = cell_properties_endo.GetLastActionPotentialDuration(90);
        double mid_APD = cell_properties_mid.GetLastActionPotentialDuration(90);

        //check that percentage increase from epi to mid and endo (roughly) matches results
        // from McIntosh et al. Card Res, 45:397-409. 200 (Figure 1 and 2)
        TS_ASSERT_DELTA((mid_APD-epi_APD)*100/epi_APD, 36.2, 2);
        TS_ASSERT_DELTA((endo_APD-epi_APD)*100/epi_APD, 8, 2);

    }

    void TestMaleckar(void) throw (Exception)
    {
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-280,
                                                                          6,
                                                                          1000,
                                                                          4.0));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);
        Maleckar2009OdeSystem atrial_ode_system(p_solver, p_stimulus);
        
        //default values
        atrial_ode_system.SetScaleFactorGks(1.0);
        atrial_ode_system.SetScaleFactorIto(1.0);
        atrial_ode_system.SetScaleFactorGkr(1.0);
        atrial_ode_system.SetScaleFactorGna(1.0);
        atrial_ode_system.SetScaleFactorAch(1e-24);
        
        // Solve and write to file for a short time
        RunOdeSolverWithIonicModel(&atrial_ode_system,
                                   25,/*end time*/
                                   "Maleckar2009",
                                   1000);
                                   
        // Check against validated data(The full AP is attached to ticket 1194)
        CheckCellModelResults("Maleckar2009");
              
        //Test the GetIIonic method against one hardcoded value.
        TS_ASSERT_DELTA(atrial_ode_system.GetIIonic(), 1.4426, 1e-3);
    }
//    Uncomment the includes for the models too
//
//    void TestOdeSolverForN98WithSimpleStimulus(void)
//    {
//        clock_t ck_start, ck_end;
//
//        // Set stimulus
//        double magnitude_stimulus = 0.0;  // uA/cm2
//        double duration_stimulus = 0.5;  // ms
//        double start_stimulus = 10.0;   // ms
//        SimpleStimulus stimulus(magnitude_stimulus,
//                                 duration_stimulus,
//                                 start_stimulus);
//
//        // Solve forward
//        HeartConfig::Instance()->SetOdeTimeStep(0.0005);
//        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
//        CML_noble_varghese_kohl_noble_1998_basic_pe_lut n98_ode_system(p_solver, p_stimulus);
//
//        std::vector<double> dY(22);
//
//        n98_ode_system.EvaluateYDerivatives(0, n98_ode_system.rGetStateVariables(), dY);
//
//        for (unsigned i = 0; i < dY.size(); ++i)
//        {
//            std::cout << dY[i] << "\n";
//        }
//
//        ck_start = clock();
//        RunOdeSolverWithIonicModel(&n98_ode_system,
//                                   150.0,
//                                   "N98RegResult");
//        ck_end = clock();
//        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
//
//        TS_ASSERT( !std::isnan( n98_ode_system.GetIIonic() ) );
//
//        // Solve backward
//        HeartConfig::Instance()->SetOdeTimeStep(0.1);
//        CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be n98_be_ode_system(p_stimulus);
//
//        ck_start = clock();
//        RunOdeSolverWithIonicModel(&n98_be_ode_system,
//                                   150.0,
//                                   "BeN98RegResult");
//        ck_end = clock();
//        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
//
//        TS_ASSERT( !std::isnan( n98_be_ode_system.GetIIonic() ) );
//
//        CompareCellModelResults("N98RegResult", "BeN98RegResult", 0.2);
//
//        std::cout << "Run times:\n\tForward: " << forward
//                  << "\n\tBackward: " << backward
//                  << std::endl;
//
//
//    }

    void TestLR1991Archiving(void) throw(Exception)
    {
        //Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename =  ArchiveLocationInfo::GetProcessUniqueFilePath("lr91.arch");

        // Save
        {
            // Set stimulus
            double magnitude_stimulus = -3;  // uA/cm2
            double duration_stimulus = 3;  // ms
            double start_stimulus = 10.0;   // ms
            boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                            duration_stimulus,
                                                                            start_stimulus));
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            double time_step = 0.01;

            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

            // Check Standard
            AbstractCardiacCell* const p_luo_rudy_cell = new LuoRudyIModel1991OdeSystem(p_solver, p_stimulus);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch <<  p_luo_rudy_cell;

            // These results are in the repository and should be replicated after the load below
//            RunOdeSolverWithIonicModel(p_luo_rudy_cell,
//                           50.0,
//                           "LRAfterArchiveValidData");

            delete p_luo_rudy_cell;
        }
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCell* p_luo_rudy_cell;
            input_arch >> p_luo_rudy_cell;

            TS_ASSERT_EQUALS( p_luo_rudy_cell->GetNumberOfStateVariables(), 8U );


            RunOdeSolverWithIonicModel(p_luo_rudy_cell,
                                       50.0,
                                       "LRAfterArchive");

            CheckCellModelResults("LRAfterArchive");

            delete p_luo_rudy_cell;
        }
     }

    void TestBackwardCellsArchiving(void) throw(Exception)
    {
        //Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename =  ArchiveLocationInfo::GetProcessUniqueFilePath("backward_cells.arch");

        // Save
        {
            // Set stimulus
            double magnitude_stimulus = -3;  // uA/cm2
            double duration_stimulus = 3;  // ms
            double start_stimulus = 10.0;   // ms
            boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                            duration_stimulus,
                                                                            start_stimulus));

            double time_step = 0.01;

            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

            AbstractCardiacCell* const p_backward_cell1 = new BackwardEulerLuoRudyIModel1991(p_stimulus);
            AbstractCardiacCell* const p_backward_cell2 = new BackwardEulerFoxModel2002Modified(p_stimulus);
            AbstractCardiacCell* const p_backward_cell3 = new BackwardEulerNobleVargheseKohlNoble1998(p_stimulus);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch <<  p_backward_cell1;
            output_arch <<  p_backward_cell2;
            output_arch <<  p_backward_cell3;

            // These results are in the repository and should be replicated after the load below
//            RunOdeSolverWithIonicModel(p_backward_cell1,
//                                       50.0,
//                                       "Backward1AfterArchiveValidData");
//
//            RunOdeSolverWithIonicModel(p_backward_cell2,
//                                       50.0,
//                                       "Backward2AfterArchiveValidData");
//
//            RunOdeSolverWithIonicModel(p_backward_cell3,
//                                       50.0,
//                                       "Backward3AfterArchiveValidData");

            delete p_backward_cell1;
            delete p_backward_cell2;
            delete p_backward_cell3;
        }
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCell* p_backward_cell1;
            AbstractCardiacCell* p_backward_cell2;
            AbstractCardiacCell* p_backward_cell3;
            input_arch >> p_backward_cell1;
            input_arch >> p_backward_cell2;
            input_arch >> p_backward_cell3;

            TS_ASSERT_EQUALS( p_backward_cell1->GetNumberOfStateVariables(), 8U );
            TS_ASSERT_EQUALS( p_backward_cell2->GetNumberOfStateVariables(), 13U );
            TS_ASSERT_EQUALS( p_backward_cell3->GetNumberOfStateVariables(), 22U );

            RunOdeSolverWithIonicModel(p_backward_cell1,
                                       50.0,
                                       "Backward1AfterArchive");

            RunOdeSolverWithIonicModel(p_backward_cell2,
                                       50.0,
                                       "Backward2AfterArchive");

            RunOdeSolverWithIonicModel(p_backward_cell3,
                                       50.0,
                                       "Backward3AfterArchive");

            CheckCellModelResults("Backward1AfterArchive");
            CheckCellModelResults("Backward2AfterArchive");
            CheckCellModelResults("Backward3AfterArchive");

            delete p_backward_cell1;
            delete p_backward_cell2;
            delete p_backward_cell3;
        }
     }

    void TestPyCMLArchiving(void) throw(Exception)
    {
        //Archive
        OutputFileHandler handler("archive", false);
        ArchiveLocationInfo::SetArchiveDirectory(handler.GetOutputDirectoryFullPath());
        std::string archive_filename =  ArchiveLocationInfo::GetProcessUniqueFilePath("noble98.arch");

        // Save
        {
            // Set stimulus
            double magnitude_stimulus = -3;  // uA/cm2
            double duration_stimulus = 3;  // ms
            double start_stimulus = 10.0;   // ms
            boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude_stimulus,
                                                                            duration_stimulus,
                                                                            start_stimulus));
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            double time_step = 0.01;

            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

            // Check Standard
            AbstractCardiacCell* const p_n98_cell = new CML_noble_varghese_kohl_noble_1998_basic(p_solver, p_stimulus);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch <<  p_n98_cell;

            // These results are in the repository and should be replicated after the load below
//            RunOdeSolverWithIonicModel(p_n98_cell,
//                                       50.0,
//                                       "N98AfterArchiveValidData");

            delete p_n98_cell;
        }
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCell* p_n98_cell;
            input_arch >> p_n98_cell;

            TS_ASSERT_EQUALS( p_n98_cell->GetNumberOfStateVariables(), 22U );

            RunOdeSolverWithIonicModel(p_n98_cell,
                                       50.0,
                                       "N98AfterArchive");

            CheckCellModelResults("N98AfterArchive");

            delete p_n98_cell;
        }
     }

private:
    void TryTestLr91WithVoltageDrop(unsigned ratio) //throw (Exception)
    {
        double end_time = 10;        // ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01/ratio, 0.01, 0.01);

        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus);
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        LuoRudyIModel1991OdeSystem lr91_ode_system(p_solver, p_zero_stimulus);
        double time=0.0;
        double start_voltage=-83.853;
        double end_voltage=-100;
        while (time<end_time)
        {
            double next_time=time + HeartConfig::Instance()->GetPdeTimeStep();
            lr91_ode_system.SetVoltage( start_voltage + (end_voltage-start_voltage)*time/end_time );
            lr91_ode_system.ComputeExceptVoltage(time, next_time);
            time=next_time;
        }
    }
};


#endif //_TESTIONICMODELS_HPP_
