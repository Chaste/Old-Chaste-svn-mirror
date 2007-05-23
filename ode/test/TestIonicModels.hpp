#ifndef _TESTIONICMODELS_HPP_
#define _TESTIONICMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "RunAndCheckIonicModels.hpp"

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

// Note: RunOdeSolverWithIonicModel(), SaveSolution(), CheckCellModelResults(), CompareCellModelResults()
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
    
    
    void TestOdeSolverForFHN61WithInitialStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude_stimulus = -80.0;   // dimensionless
        double duration_stimulus = 0.5 ;  // ms
        double start_stimulus = 0.0;   // ms
        InitialStimulus stimulus(magnitude_stimulus,
                                 duration_stimulus,
                                 start_stimulus);
                                 
        EulerIvpOdeSolver solver;
        double time_step = 0.01;
        FitzHughNagumo1961OdeSystem fhn61_ode_system(&solver, time_step, &stimulus);
        
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
        double duration  = 2.0  ;  // ms
        double when = 50.0; // ms
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds
        
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, time_step, &stimulus);
        
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
        double duration  = 2.0  ;  // ms
        double start = 50.0; // ms
        double frequency = 1.0/500; // ms^-1
        RegularStimulus stimulus(magnitude, duration, frequency, start);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds
        
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, time_step, &stimulus);
        
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
        double duration  = 2.0  ;  // ms
        double when = 50.0; // ms
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds
        
        // Solve using backward euler
        BackwardEulerLuoRudyIModel1991 lr91_backward_euler(time_step, &stimulus);
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
                                   "Lr91DelayedSt= i_stim im");
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        
        // Compare results
        CompareCellModelResults("Lr91DelayedStim", "Lr91BackwardEuler", 0.01);
        
        // Try with larger timestep and coarser tolerance.
        // We can't use a larger time step than 0.01 for forward Euler - the gating
        // variables go out of range.
        BackwardEulerLuoRudyIModel1991 lr91_backward_euler2(time_step*50, &stimulus);
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
        LuoRudyIModel1991OdeSystem lr91(&solver, 0.01, &stimulus);
        BackwardEulerLuoRudyIModel1991 backward_lr91(0.01, &stimulus);
        // calc IIonic using initial conditions
        TS_ASSERT_DELTA(lr91.GetIIonic(), backward_lr91.GetIIonic(), 1e-12);
    }
    
    void TestOdeSolverForFox2002WithRegularStimulus(void) throw (Exception)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -80.0;
        double duration  = 1.0  ;  // ms
        double start = 50.0; // ms
        double frequency = 1.0/500; // ms^-1
        RegularStimulus stimulus(magnitude, duration, frequency, start);
        
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
        
        CompareCellModelResults("FoxRegularStim", "BackwardFoxRegularStim", 0.15);
        
        std::cout << "Run times:\n\tForward: " << forward
                  << "\n\tBackward: " << backward
                  << std::endl;
        
    }
    
    void TestLr91WithVoltageDropVariousTimeStepRatios()
    {
        TS_ASSERT_THROWS_ANYTHING(TryTestLr91WithVoltageDrop(1));
        TS_ASSERT_THROWS_ANYTHING(TryTestLr91WithVoltageDrop(2));
        TS_ASSERT_THROWS_ANYTHING(TryTestLr91WithVoltageDrop(3));
        TS_ASSERT_THROWS_NOTHING(TryTestLr91WithVoltageDrop(4));
           
    }
    
private:
    void TryTestLr91WithVoltageDrop(unsigned ratio) throw (Exception)
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
