#ifndef _TESTIONICMODELS_HPP_
#define _TESTIONICMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

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

const double TOLERANCE = 1e-2; // Not used at present


class TestIonicModels : public CxxTest::TestSuite
{
public:

    void RunOdeSolverWithIonicModel(AbstractCardiacCell *pOdeSystem,
                                    double endTime,
                                    const char *pFilename)
    {
        double start_time = 0.0;
        
        // Store the current system state
        std::vector<double> state_variables_ref = pOdeSystem->rGetStateVariables();
        std::vector<double> state_variables_copy = state_variables_ref;
        
        // Test ComputeExceptVoltage
        double v_init = pOdeSystem->GetVoltage();
        OdeSolution solution = pOdeSystem->ComputeExceptVoltage(start_time, endTime);
        double v_end = pOdeSystem->GetVoltage();
        TS_ASSERT_DELTA(v_init, v_end, 1e-6);
        
        // Test SetVoltage
        pOdeSystem->SetVoltage(1e6);
        TS_ASSERT_DELTA(pOdeSystem->GetVoltage(), 1e6, 1e-6);
        
        // Reset the system
        pOdeSystem->SetStateVariables(state_variables_copy);
        
        /*
         * Solve 
         */
        solution = pOdeSystem->Compute(start_time, endTime);
        
        /*
         * Write data to a file using ColumnDataWriter
         */
        int step_per_row = 100;
        ColumnDataWriter writer("TestIonicModels",pFilename,false);
        int time_var_id = writer.DefineUnlimitedDimension("Time","ms");
        
        std::vector<int> var_ids;
        for (unsigned i=0; i<pOdeSystem->rGetVariableNames().size(); i++)
        {
            var_ids.push_back(writer.DefineVariable(pOdeSystem->rGetVariableNames()[i],
                                                    pOdeSystem->rGetVariableUnits()[i]));
        }
        writer.EndDefineMode();
        
        for (unsigned i = 0; i < solution.rGetSolutions().size(); i+=step_per_row)
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
        /*
         * Compare 2 sets of results, e.g. from 2 different solvers for the same model.
         * Initially we assume the time series are the same; this will change.
         */
        
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
        
        TS_ASSERT_EQUALS(times1.size(), times2.size());
        for (unsigned i=0; i<times1.size(); i++)
        {
            TS_ASSERT_DELTA(times1[i], times2[i], 1e-12);
            // adjust tol to data
            TS_ASSERT_DELTA(voltages1[i], voltages2[i], tolerance);
            TS_ASSERT_DELTA(calcium1[i],  calcium2[i],  tolerance/1000);
            TS_ASSERT_DELTA(h1[i],        h2[i],        tolerance/10);
        }
    }
    
    void testOdeSolverForHH52WithInitialStimulus(void)
    {
        /*
         * Set stimulus
         */
        double magnitude_stimulus = 20.0;  // uA/cm2
        double duration_stimulus = 0.5;  // ms
        double start_stimulus = 10.0;   // ms
        InitialStimulus stimulus(magnitude_stimulus,
                                 duration_stimulus,
                                 start_stimulus);
        EulerIvpOdeSolver solver;
        double time_step = 0.01;
        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(&solver, time_step, &stimulus);
        
        /*
         * Solve and write to file
         */
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
    
    
    void testOdeSolverForFHN61WithInitialStimulus(void)
    {
        /*
         * Set stimulus
         */
        double magnitude_stimulus = -80.0;   // dimensionless
        double duration_stimulus = 0.5 ;  // ms
        double start_stimulus = 0.0;   // ms
        InitialStimulus stimulus(magnitude_stimulus,
                                 duration_stimulus,
                                 start_stimulus);
                                 
        EulerIvpOdeSolver solver;
        double time_step = 0.01;
        FitzHughNagumo1961OdeSystem fhn61_ode_system(&solver, time_step, &stimulus);
        
        /*
         * Solve and write to file
         */
        RunOdeSolverWithIonicModel(&fhn61_ode_system,
                                   500.0,
                                   "FHN61RegResult");
                                   
        CheckCellModelResults("FHN61RegResult");
        
        // test GetIionic ('fake' ionic current) (the GetIionic method was first
        // manually tested by changing the EvaluateYDerivatives() code to call it,
        // this verified that GetIionic has no errors, therefore we can test here
        // against a hardcoded result
        TS_ASSERT_DELTA( fhn61_ode_system.GetIIonic(), -0.0058, 1e-3);
    }
    
    
    void testOdeSolverForLR91WithDelayedInitialStimulus(void)
    {
        /*
         * Set stimulus
         */
        double magnitude = -25.5;
        double duration  = 2.0  ;  // ms
        double when = 50.0; // ms
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds
        
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, time_step, &stimulus);
        
        /*
         * Solve and write to file
         */
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
    
    void testOdeSolverForLR91WithRegularStimulus(void) throw (Exception)
    {
        /*
         * Set stimulus
         */
        double magnitude = -25.5;
        double duration  = 2.0  ;  // ms
        double start = 50.0; // ms
        double frequency = 1.0/500; // ms^-1
        RegularStimulus stimulus(magnitude, duration, frequency, start);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds
        
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, time_step, &stimulus);
        
        /*
         * Solve and write to file
         */
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91RegularStim");
                                   
        CheckCellModelResults("Lr91RegularStim");
        
    }
    
    void testBackwardEulerLr91WithDelayedInitialStimulus(void) throw (Exception)
    {
        /*
         * Set stimulus
         */
        double magnitude = -25.5;
        double duration  = 2.0  ;  // ms
        double when = 50.0; // ms
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds
        
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, time_step, &stimulus);
        
        /*
         * Solve and write to file
         */
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91DelayedStim");
                                   
        // Now use the backward euler solver
        BackwardEulerLuoRudyIModel1991 lr91_backward_euler(time_step, &stimulus);
        RunOdeSolverWithIonicModel(&lr91_backward_euler,
                                   end_time,
                                   "Lr91BackwardEuler");
        
        // Compare results
        CompareCellModelResults("Lr91DelayedStim", "Lr91BackwardEuler",0.25);
    }
};


#endif //_TESTIONICMODELS_HPP_
