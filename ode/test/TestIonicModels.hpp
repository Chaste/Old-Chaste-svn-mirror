#ifndef _TESTIONICMODELS_HPP_
#define _TESTIONICMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AdamsBashforthIvpOdeSolver.hpp"

#include "OdeSolution.hpp"

#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"

#include "AbstractOdeSystem.hpp"
#include "AbstractCardiacCell.hpp"
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "Noble98ForwardEulerFromCellml.hpp"
#include "Noble98BackwardEulerFromCellml.hpp"


class TestIonicModels : public CxxTest::TestSuite
{
public:

    void RunOdeSolverWithIonicModel(AbstractCardiacCell *pOdeSystem,
                                    double endTime,
                                    std::string filename,
                                    int stepPerRow=100,
                                    bool doComputeExceptVoltage=true)
    {
        double start_time = 0.0;
        
        if (doComputeExceptVoltage)
        {
            // Store the current system state
            std::vector<double> state_variables_ref = pOdeSystem->rGetStateVariables();
            std::vector<double> state_variables_copy = state_variables_ref;
            
            // Test ComputeExceptVoltage
            double v_init = pOdeSystem->GetVoltage();
            OdeSolution solution = pOdeSystem->ComputeExceptVoltage(start_time, endTime);
            double v_end = pOdeSystem->GetVoltage();
            TS_ASSERT_DELTA(v_init, v_end, 1e-6);
            
            // Save results for comparison
            // This appears to be the only use of the return value of ComputeExceptVoltage
            SaveSolution(filename + "_ExceptVoltage", pOdeSystem, solution, stepPerRow);
            
            // Test SetVoltage
            pOdeSystem->SetVoltage(1e6);
            TS_ASSERT_DELTA(pOdeSystem->GetVoltage(), 1e6, 1e-6);
            
            // Reset the system
            pOdeSystem->SetStateVariables(state_variables_copy);
        }
        
        // Solve
        OdeSolution solution = pOdeSystem->Compute(start_time, endTime);
        
        // Write data to a file using ColumnDataWriter
        SaveSolution(filename, pOdeSystem, solution, stepPerRow);
    }
    
    void SaveSolution(std::string baseResultsFilename, AbstractCardiacCell *pOdeSystem,
                      OdeSolution& solution, int stepPerRow)
    {
        // Write data to a file using ColumnDataWriter
        
        ColumnDataWriter writer("TestIonicModels",baseResultsFilename,false);
        int time_var_id = writer.DefineUnlimitedDimension("Time","ms");
        
        std::vector<int> var_ids;
        for (unsigned i=0; i<pOdeSystem->rGetVariableNames().size(); i++)
        {
            var_ids.push_back(writer.DefineVariable(pOdeSystem->rGetVariableNames()[i],
                                                    pOdeSystem->rGetVariableUnits()[i]));
        }
        writer.EndDefineMode();
        
        for (unsigned i = 0; i < solution.rGetSolutions().size(); i+=stepPerRow)
        {
            writer.PutVariable(time_var_id, solution.rGetTimes()[i]);
            for (unsigned j=0; j<var_ids.size(); j++)
            {
                writer.PutVariable(var_ids[j], solution.rGetSolutions()[i][j]);
            }
            writer.AdvanceAlongUnlimitedDimension();
        }
        writer.Close();
    }
    
    void CheckCellModelResults(std::string baseResultsFilename)
    {
        /*
         * Check the cell model against a previous version
         * and another source e.g. Alan's COR
         */
        
        // read data entries for the new file and compare to valid data from
        // other source
        ColumnDataReader data_reader("TestIonicModels", baseResultsFilename);
        std::vector<double> times = data_reader.GetValues("Time");
        std::vector<double> voltages = data_reader.GetValues("V");
        ColumnDataReader valid_reader("ode/test/data", baseResultsFilename+"ValidData",
                                      false);
        std::vector<double> valid_times = valid_reader.GetValues("Time");
        std::vector<double> valid_voltages = valid_reader.GetValues("V");
        
        TS_ASSERT_EQUALS(times.size(), valid_times.size());
        for (unsigned i=0; i<valid_times.size(); i++)
        {
            TS_ASSERT_DELTA(times[i], valid_times[i], 1e-12);
            // adjust tol to data
            TS_ASSERT_DELTA(voltages[i], valid_voltages[i], 1e-6);
        }
    }
    
    void CompareCellModelResults(std::string baseResultsFilename1, std::string baseResultsFilename2, double tolerance)
    {
        // Compare 2 sets of results, e.g. from 2 different solvers for the same model.
        // Initially we assume the time series are the same; this will change.
        // If the time series differ, the finer resolution must be given first.
        std::cout << "Comparing " << baseResultsFilename1 << " with "
        << baseResultsFilename2 << std::endl;
        
        ColumnDataReader data_reader1("TestIonicModels", baseResultsFilename1);
        std::vector<double> times1 = data_reader1.GetValues("Time");
        std::vector<double> voltages1 = data_reader1.GetValues("V");
        std::vector<double> calcium1 = data_reader1.GetValues("CaI");
        std::vector<double> h1 = data_reader1.GetValues("h");
        
        ColumnDataReader data_reader2("TestIonicModels", baseResultsFilename2);
        std::vector<double> times2 = data_reader2.GetValues("Time");
        std::vector<double> voltages2 = data_reader2.GetValues("V");
        std::vector<double> calcium2 = data_reader2.GetValues("CaI");
        std::vector<double> h2 = data_reader2.GetValues("h");
        
        TS_ASSERT(times1.size() >= times2.size());
        double last_v = voltages2[0];
        double tol = tolerance;
        for (unsigned i=0, j=0; i<times2.size(); i++)
        {
            // Find corresponding time index
            while (j<times1.size() && times1[j] < times2[i] - 1e-12)
            {
                j++;
            }
            
            // Set tolerance higher in upstroke
            if (fabs(voltages2[i] - last_v) > 0.05)
            {
                tol = tolerance * 25;
            }
            else
            {
                tol = tolerance;
            }
            last_v = voltages2[i];
            
            TS_ASSERT_DELTA(times1[j], times2[i], 1e-12);
            // adjust tol to data
            TS_ASSERT_DELTA(voltages1[j], voltages2[i], tol);
            TS_ASSERT_DELTA(calcium1[j],  calcium2[i],  tol/1000);
            TS_ASSERT_DELTA(h1[j],        h2[i],        tol/10);
        }
    }
    
    void TestOdeSolverForHH52WithInitialStimulus(void)
    {
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
        RunOdeSolverWithIonicModel(&hh52_ode_system,
                                   150.0,
                                   "HH52RegResult");
                                   
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
        RunOdeSolverWithIonicModel(&fhn61_ode_system,
                                   500.0,
                                   "FHN61RegResult");
                                   
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
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91DelayedStim");
                                   
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
                                   "Lr91DelayedStim");
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
    
    void TestNoble98WithDelayedInitialStimulus(void) throw (Exception)
    {
        clock_t ck_start, ck_end;
        // Set stimulus
        double magnitude = -3;    // nA
        double duration  = 2.0/1000;   // s
        double when = 60.0/1000;       // s
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 0.2;   // seconds
        double time_step = 2e-6; // seconds
        
        // Solve using forward euler.
        EulerIvpOdeSolver solver;
        CML_noble_model_1998 n98_forward_euler(&solver, time_step, &stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_forward_euler,
                                   end_time,
                                   "N98ForwardEuler", 300, false);
        ck_end = clock();
        double forward1 = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        
        // Solve using backward euler.
        // I think a small timestep is needed because of the forward Euler step for V,
        // which means the numerical method isn't stable for large dt.
        // This is annoying.
        CML_noble_model_1998_BE n98_backward_euler(time_step, &stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&n98_backward_euler,
                                   end_time,
                                   "N98BackwardEuler", 300, false);
        ck_end = clock();
        double backward1 = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        
        CheckCellModelResults("N98BackwardEuler");
        CompareCellModelResults("N98ForwardEuler", "N98BackwardEuler", 0.015);
        
        // Check GetIIonic
        double i_ion = n98_backward_euler.GetIIonic();
        TS_ASSERT_DELTA(i_ion, 0.0228, 0.0001);
        TS_ASSERT_DELTA(i_ion, n98_forward_euler.GetIIonic(), 0.0001);
        
//        // Try with smaller timestep.
//        CML_noble_model_1998_BE n98_backward_euler2(time_step/2, &stimulus);
//        ck_start = clock();
//        RunOdeSolverWithIonicModel(&n98_backward_euler2,
//                                   end_time,
//                                   "N98BackwardEuler2", 600, false);
//        ck_end = clock();
//        double backward2 = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
//        CompareCellModelResults("N98BackwardEuler", "N98BackwardEuler2", 0.015);
//
        std::cout << "Run times:\n\tForward: " << forward1
        << "\n\tBackward: " << backward1
//            << "\n\tBackward (half dt): " << backward2
        << std::endl;
    }
};


#endif //_TESTIONICMODELS_HPP_
