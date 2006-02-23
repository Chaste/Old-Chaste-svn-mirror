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
                                            timeStep, pOdeSystem->mInitialConditions);
        
        /*
         * Write data to a file using ColumnDataWriter
         */                                                           
        ColumnDataWriter writer("testoutput",pFilename);
        int time_var_id = writer.DefineUnlimitedDimension("Time","ms");
        
        std::vector<int> var_ids;
        for (int i=0; i<pOdeSystem->mVariableNames.size(); i++)
        {
            var_ids.push_back(writer.DefineVariable(pOdeSystem->mVariableNames[i],
                                                    pOdeSystem->mVariableUnits[i]));
        }
        writer.EndDefineMode();
                
        for (int i = 0; i < solution.mSolutions.size(); i+=step_per_row) 
        {
            if (i!=0)
            {
                writer.AdvanceAlongUnlimitedDimension();
            }
            writer.PutVariable(time_var_id, solution.mTime[i]);
            for (int j=0; j<var_ids.size(); j++)
            {
                writer.PutVariable(var_ids[j], solution.mSolutions[i][j]);
            }
        }        
        writer.Close();        
    }
   
        
    void testOdeSolverForHH52WithRegularStimulus(void)
    {
        /*
         * Set stimulus
         */   
        double magnitude_of_stimulus = -20.0;  
        double duration_of_stimulus  = 0.5 ;  // ms                     
        double frequency = 1.0/50.0;
        double when = 40.0;                                      
        RegularStimulus stimulus(magnitude_of_stimulus,
                                 duration_of_stimulus,
                                 frequency,
                                 when);

        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(&stimulus);
        
        /*
         * Solve and write to file
         */
        runOdeSolverWithIonicModel(&hh52_ode_system,
                                   150.0,
                                   0.01,
                                   "HH52RegResult");
    }


    void testOdeSolverForFHN61WithRegularStimulus(void)
    {
        /*
         * Set stimulus
         */             
        double magnitude = 1.0;  
        double duration  = 0.5 ;  // ms                     
        double frequency = 1.0/500.0;
        double start_stimulus = 40.0;                 
        RegularStimulus stimulus(magnitude, duration, frequency, start_stimulus); 

        FitzHughNagumo1961OdeSystem fhn61_ode_system(&stimulus);
        
        /*
         * Solve and write to file
         */
        runOdeSolverWithIonicModel(&fhn61_ode_system,
                                   500.0,
                                   0.01,
                                   "FHN61RegResult");
    }
    
        
    void testOdeSolverForLR91WithDelayedInitialStimulus(void)
    {
        /*
         * Set stimulus
         */
        double magnitude = -80.0;  
        double duration  = 0.5 ;  // ms                     
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
                                   "NewDelayedStimLR91");
        
        //read in good data file and compare line by line
        std::ifstream testfile("testoutput/NewDelayedStimLR91.dat",std::ios::in);
        std::ifstream goodfile("ode/test/data/Lr91DelayedStimGood.dat",std::ios::in);

        std::string teststring;
        std::string goodstring;
        
        while(getline(testfile, teststring))
        {
              getline(goodfile,goodstring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();                       

        
        // read data entries for the new file and compare to valid data from 
        // ????????????        
        ColumnDataReader data_reader("testoutput","NewDelayedStimLR91");
        std::vector<double> times = data_reader.GetValues("Time");
        std::vector<double> voltages = data_reader.GetValues("V");

        ColumnDataReader valid_reader("ode/test/data/","Lr91DelayedStimValidData");
        std::vector<double> valid_times = valid_reader.GetValues("Time");
        std::vector<double> valid_voltages = valid_reader.GetValues("V");
        
        for(int i=0; i<valid_times.size(); i++)
        {
            TS_ASSERT_DELTA(times[i], valid_times[i], 1e-6);
            // adjust tol to data
            TS_ASSERT_DELTA(voltages[i], valid_voltages[i], 1e-6);
        }
    }   
};


#endif //_TESTIONICMODELS_HPP_
