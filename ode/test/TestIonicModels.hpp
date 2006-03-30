#ifndef _TESTIONICMODELS_HPP_
#define _TESTIONICMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AdamsBashforthIvpOdeSolver.hpp"

#include "OdeSolution.hpp"

#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"

#include "AbstractOdeSystem.hpp"
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

const double TOLERANCE = 1e-2; // Not used at present


class TestIonicModels : public CxxTest::TestSuite
{
public:

    void runOdeSolverWithIonicModel(AbstractOdeSystem *pOdeSystem,
                                    double endTime,
                                    double timeStep,
                                    const char *pFilename)
    {
        EulerIvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double start_time = 0.0;
        int step_per_row = 100;             
                
        OdeSolution solution = solver.Solve(pOdeSystem, start_time, endTime,
                                            timeStep, pOdeSystem->GetInitialConditions());
        
        /*
         * Write data to a file using ColumnDataWriter
         */                                                           
        ColumnDataWriter writer("testoutput",pFilename);
        int time_var_id = writer.DefineUnlimitedDimension("Time","ms");
        
        std::vector<int> var_ids;
        for (unsigned i=0; i<pOdeSystem->mVariableNames.size(); i++)
        {
            var_ids.push_back(writer.DefineVariable(pOdeSystem->mVariableNames[i],
                                                    pOdeSystem->mVariableUnits[i]));
        }
        writer.EndDefineMode();
                
        for (unsigned i = 0; i < solution.mSolutions.size(); i+=step_per_row) 
        {
            if (i!=0)
            {
                writer.AdvanceAlongUnlimitedDimension();
            }
            writer.PutVariable(time_var_id, solution.mTime[i]);
            for (unsigned j=0; j<var_ids.size(); j++)
            {
                writer.PutVariable(var_ids[j], solution.mSolutions[i][j]);
            }
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
        ColumnDataReader data_reader("testoutput",baseResultsFilename);
        std::vector<double> times = data_reader.GetValues("Time");
        std::vector<double> voltages = data_reader.GetValues("V");
        ColumnDataReader valid_reader("ode/test/data",baseResultsFilename+"ValidData");
        std::vector<double> valid_times = valid_reader.GetValues("Time");
        std::vector<double> valid_voltages = valid_reader.GetValues("V");
       
        for(unsigned i=0; i<valid_times.size(); i++)
        {
            TS_ASSERT_DELTA(times[i], valid_times[i], 1e-6);
            // adjust tol to data
            TS_ASSERT_DELTA(voltages[i], valid_voltages[i], 1e-6);
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

        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(&stimulus);
        
        /*
         * Solve and write to file
         */
        runOdeSolverWithIonicModel(&hh52_ode_system,
                                   150.0,
                                   0.01,
                                   "HH52RegResult");
                                   
        CheckCellModelResults("HH52RegResult");
    }


    void testOdeSolverForFHN61WithInitialStimulus(void)
    {
        /*
         * Set stimulus
         */             
        double magnitude_stimulus = -80.0;   // dimensionless
        double duration_stimulus = 0.499999 ;  // ms     //\todo Alan to regenerate                
        double start_stimulus = 0.0;   // ms
        InitialStimulus stimulus(magnitude_stimulus,
                                 duration_stimulus,
                                 start_stimulus); 

        FitzHughNagumo1961OdeSystem fhn61_ode_system(&stimulus);
        
        /*
         * Solve and write to file
         */
        runOdeSolverWithIonicModel(&fhn61_ode_system,
                                   500.0,
                                   0.01,
                                   "FHN61RegResult");
                                   
         CheckCellModelResults("FHN61RegResult");

    }
    
        
    void testOdeSolverForLR91WithDelayedInitialStimulus(void)
    {
        /*
         * Set stimulus
         */
        double magnitude = -80.0;  
        double duration  = 0.49999  ;  // ms                           
        double when = 100.0; // ms
        InitialStimulus stimulus(magnitude, duration, when); 
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.01;  //1e-5 seconds in milliseconds           
        
        LuoRudyIModel1991OdeSystem lr91_ode_system(&stimulus);
        
        /*
         * Solve and write to file
         */
        runOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   time_step,
                                   "Lr91DelayedStim");

         CheckCellModelResults("Lr91DelayedStim");
    }   
};


#endif //_TESTIONICMODELS_HPP_
