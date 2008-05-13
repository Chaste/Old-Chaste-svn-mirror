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


#ifndef _TESTIONICMODELS_HPP_
#define _TESTIONICMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "RunAndCheckIonicModels.hpp"
#include "Exception.hpp"

#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BackwardEulerLuoRudyIModel1991.hpp"

#include "FoxModel2002Modified.hpp"
#include "BackwardEulerFoxModel2002Modified.hpp"

#include "FaberRudy2000Version3.cpp"
#include "FaberRudy2000Version3Optimised.hpp"

// Note: RunOdeSolverWithIonicModel(), CheckCellModelResults(), CompareCellModelResults()
// are defined in RunAndCheckIonicModels.hpp

class TestIonicModels : public CxxTest::TestSuite
{
public:
    void TestOdeSolverForHH52WithInitialStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = 20.0;  // uA/cm2
        double duration_stimulus = 0.5;  // ms
        double start_stimulus = 10.0;   // ms
        InitialStimulus stimulus(magnitude_stimulus,
                                 duration_stimulus,
                                 start_stimulus);
        EulerIvpOdeSolver solver;
        double time_step = 0.01;
        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(&solver, time_step, &stimulus);
        
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
    
    
    void TestOdeSolverForFHN61WithInitialStimulus(void) throw (Exception)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -80.0;   // dimensionless
        double duration_stimulus = 0.5;  // ms
        double start_stimulus = 0.0;   // ms
        InitialStimulus stimulus(magnitude_stimulus,
                                 duration_stimulus,
                                 start_stimulus);
                                 
        EulerIvpOdeSolver solver;
        double time_step = 0.01;
        FitzHughNagumo1961OdeSystem fhn61_ode_system(&solver, time_step, &stimulus);

        // fhn has no [Ca_i]
        TS_ASSERT_THROWS_ANYTHING(fhn61_ode_system.GetIntracellularCalciumConcentration());
        
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
        InitialStimulus another_stimulus(-200,1.0, 0.0);
        InitialStimulus intra_stimulus(-100,1.0, 0.0);
        InitialStimulus extra_stimulus(-50, 1.0, 0.0);
        FitzHughNagumo1961OdeSystem another_fhn61_ode_system(&solver, time_step, &stimulus);
        
        another_fhn61_ode_system.SetStimulusFunction(&another_stimulus);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetStimulus(0.5), -200, 1e-12);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetIntracellularStimulus(0.5), -200, 1e-12);
        
        another_fhn61_ode_system.SetIntracellularStimulusFunction(&intra_stimulus);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetStimulus(0.5), -100, 1e-12);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetIntracellularStimulus(0.5), -100, 1e-12);
        
        another_fhn61_ode_system.SetExtracellularStimulusFunction(&extra_stimulus);
        TS_ASSERT_DELTA(another_fhn61_ode_system.GetExtracellularStimulus(0.5), -50, 1e-12);
    }
    
    
    void TestOdeSolverForLR91WithDelayedInitialStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds
        
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, time_step, &stimulus);
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
        RegularStimulus stimulus(magnitude, duration, period, start);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds
        
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, time_step, &stimulus);
        
        // cover get intracellular calcium
        TS_ASSERT_DELTA(lr91_ode_system.GetIntracellularCalciumConcentration(), 0.0002, 1e-5)
        
        // Solve and write to file
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91RegularStim");
                                   
        CheckCellModelResults("Lr91RegularStim");
    }
    
    void TestBackwardEulerLr91WithDelayedInitialStimulus(void) throw (Exception)
    {
        clock_t ck_start, ck_end;
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds
        
        // Solve using backward euler
        BackwardEulerLuoRudyIModel1991 lr91_backward_euler(time_step, &stimulus);
        
        // cover get intracellular calcium
        TS_ASSERT_DELTA(lr91_backward_euler.GetIntracellularCalciumConcentration(), 0.0002, 1e-5)
        
        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_backward_euler,
                                   end_time,
                                   "Lr91BackwardEuler");
        ck_end = clock();
        double backward1 = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        
        // Solve using forward Euler
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, time_step, &stimulus);
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
        EulerIvpOdeSolver* p_solver = new EulerIvpOdeSolver();
        BackwardEulerLuoRudyIModel1991 lr91_backward_euler2(p_solver ,time_step*50, &stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_backward_euler2,
                                   end_time,
                                   "Lr91BackwardEuler2");
        delete p_solver;
        ck_end = clock();
        double backward2 = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        CompareCellModelResults("Lr91DelayedStim", "Lr91BackwardEuler2", 0.25);
        
        std::cout << "Run times:\n\tForward: " << forward << "\n\tBackward: "
                  << backward1 << "\n\tBackward (long dt): " << backward2 << std::endl;
        
        
        // cover and check GetIIonic() match for normal and backward euler lr91
        LuoRudyIModel1991OdeSystem lr91(&solver, 0.01, &stimulus);
        BackwardEulerLuoRudyIModel1991 backward_lr91(0.01, &stimulus);
        // calc IIonic using initial conditions
        TS_ASSERT_DELTA(lr91.GetIIonic(), backward_lr91.GetIIonic(), 1e-12);

        // cover alternative constructor
        EulerIvpOdeSolver solver2;
        BackwardEulerLuoRudyIModel1991 lr91_backward_euler3(&solver2, time_step, &stimulus);
    }
    
    void TestOdeSolverForFR2000WithDelayedInitialStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0;  // ms
        double when = 10.0; // ms
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 1000.0; //ms
        double time_step = 0.007;
        
        EulerIvpOdeSolver solver;
        FaberRudy2000Version3Optimised fr2000_ode_system_opt(&solver, time_step, &stimulus);
        FaberRudy2000Version3 fr2000_ode_system(&solver, time_step, &stimulus);
        
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
    }
    
    
    void TestOdeSolverForFR2000WithVariablePotassiumCurrents(void)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0;  // ms
        double when = 0.0; // ms
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 1000.0; //ms
        double time_step = 0.007;
        
        EulerIvpOdeSolver solver;
        FaberRudy2000Version3 fr2000_ode_system_endo(&solver, time_step, &stimulus);
        fr2000_ode_system_endo.SetScaleFactorGks(0.462);
        fr2000_ode_system_endo.SetScaleFactorIto(0.0);
        
        // Solve and write to file
        RunOdeSolverWithIonicModel(&fr2000_ode_system_endo,
                                   end_time,
                                   "FR2000Endo",
                                   500, false);
        
        CheckCellModelResults("FR2000Endo");
        
        FaberRudy2000Version3 fr2000_ode_system_mid(&solver, time_step, &stimulus);
        fr2000_ode_system_mid.SetScaleFactorGks(1.154);
        fr2000_ode_system_mid.SetScaleFactorIto(0.85);
        
        // Solve and write to file
        RunOdeSolverWithIonicModel(&fr2000_ode_system_mid,
                                   end_time,
                                   "FR2000Mid",
                                   500, false);
        
        CheckCellModelResults("FR2000Mid");
        
        FaberRudy2000Version3 fr2000_ode_system_epi(&solver, time_step, &stimulus);
        fr2000_ode_system_epi.SetScaleFactorGks(1.154);
        fr2000_ode_system_epi.SetScaleFactorIto(1.0);
        
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
        RegularStimulus stimulus(magnitude, duration, period, start);
        
        double end_time = 200.0;  // milliseconds
        double time_step = 0.002; //2e-6 seconds in milliseconds
                            // 0.005 leads to NaNs.
        
        EulerIvpOdeSolver solver;
        FoxModel2002Modified fox_ode_system(&solver, time_step, &stimulus);
        BackwardEulerFoxModel2002Modified backward_system(time_step*5, &stimulus);
        
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
        TS_ASSERT_THROWS_ANYTHING(TryTestLr91WithVoltageDrop(1))
        TS_ASSERT_THROWS_ANYTHING(TryTestLr91WithVoltageDrop(2))
        TS_ASSERT_THROWS_ANYTHING(TryTestLr91WithVoltageDrop(3))
        TS_ASSERT_THROWS_NOTHING(TryTestLr91WithVoltageDrop(4));
    }
    
private:
    void TryTestLr91WithVoltageDrop(unsigned ratio) //throw (Exception)
    {
        double pde_time_step = 0.01;  // ms (not used, but here to replicate TestMonodomainHeart)
        double ode_time_step = pde_time_step/ratio; // ms
        double end_time = 10;        // ms
        InitialStimulus zero_stimulus(0,0,0);
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, ode_time_step, &zero_stimulus);
        double time=0.0;
        double start_voltage=-83.853;
        double end_voltage=-100;
        while (time<end_time)
        {   
            double next_time=time+pde_time_step;
            lr91_ode_system.SetVoltage( start_voltage + (end_voltage-start_voltage)*time/end_time );
            lr91_ode_system.ComputeExceptVoltage(time, next_time);
            time=next_time;
        }
    }
};


#endif //_TESTIONICMODELS_HPP_
