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


class TestOdeSolverForHH52 : public CxxTest::TestSuite
{
public:

    void runOdeSolverWithIonicModel(AbstractOdeSystem *pOdeSystem,
                                    double endTime,
                                    double timeStep,
                                    std::vector<double> &rInitialConditions,
                                    const char *pFilename,
                                    std::vector<std::string> &rVariableNames,
                                    std::vector<std::string> &rVariableUnits)
    {
        // Consistency checks
        assert(rInitialConditions.size() == rVariableNames.size());
        assert(rInitialConditions.size() == rVariableUnits.size());
        
        /*
         * Choose an ode solver
         */
        EulerIvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double start_time = 0.0;
        int step_per_row = 100;             
                
        OdeSolution solution = solver.Solve(pOdeSystem, start_time, endTime,
                                            timeStep, rInitialConditions);
        
        /*
         * Write data to a file using ColumnDataWriter
         */                                                           
        ColumnDataWriter writer("testoutput",pFilename);
        int time_var_id = writer.DefineUnlimitedDimension("Time","ms");
        
        std::vector<int> var_ids;
        for (int i=0; i<rVariableNames.size(); i++)
        {
            var_ids.push_back(writer.DefineVariable(rVariableNames[i],
                                                    rVariableUnits[i]));
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
   
    /* Test Ode Solver for HH52
     *
     * The aim of this test is to simulate a single cell.
     * We run the ode solver with the cell model for a period of
     * time.  
     * 
     */
    void runOdeSolverForHH52(AbstractStimulusFunction *pStimulus,
                             const char *pFilename,
                             double endTime,
                             double timeStep)
    {
        /*
         * Instantiate the ionic model: need to pass stimulus function
         */        
        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(pStimulus);
        
        /*
         * Create vectors of variable names & units
         */
        std::vector<std::string> variable_names;
        std::vector<std::string> variable_units;
        std::vector<double> initial_conditions;
        
        variable_names.push_back("V");
        variable_units.push_back("mV");
        initial_conditions.push_back(0.0);
        
        variable_names.push_back("n");
        variable_units.push_back("");
        initial_conditions.push_back(0.31768);
        
        variable_names.push_back("h");
        variable_units.push_back("");
        initial_conditions.push_back(0.59612);
        
        variable_names.push_back("m");
        variable_units.push_back("");
        initial_conditions.push_back(0.05293);
        
        /*
         * Solve and write to file
         */
        runOdeSolverWithIonicModel(&hh52_ode_system,
                                   endTime,
                                   timeStep,
                                   initial_conditions,
                                   pFilename,
                                   variable_names,
                                   variable_units);
    }
    
    void testOdeSolverForHH52WithInitialStimulus(void)
    {
        /*
         * Set magnitude of stimulus
         * 
         */
  
        double magnitude_of_stimulus = -20.0;  
        double duration_of_stimulus  = 0.5 ;  // ms                     
        
        /*
         * Choose function for stimulus
         */             
        InitialStimulus stimulus(magnitude_of_stimulus, duration_of_stimulus);


        runOdeSolverForHH52(&stimulus, "HH52Result.dat", 5.0, 0.0001);
    }	
    
    void testOdeSolverForHH52WithRegularStimulus(void)
    {
        /*
         * Set magnitude of stimulus
         * 
         */
           
        double magnitude_of_stimulus = -20.0;  
        double duration_of_stimulus  = 0.5 ;  // ms                     
        double frequency = 1.0/50.0;
        double start_stimulus = 40.0;              
        
        /*
         * Choose function for stimulus
         */                    
        RegularStimulus stimulus(magnitude_of_stimulus,
                                 duration_of_stimulus,
                                 frequency,
                                 start_stimulus);

        /*
         * Solve 
         */            
        runOdeSolverForHH52(&stimulus, "HH52RegResult.dat", 150.0, 0.01);
             
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    void runOdeSolverForFHN61(AbstractStimulusFunction *pStimulus,
                              const char *pFilename,
                              double endTime,
                              double timeStep)
    {
        /*
         * Instantiate the FHN model: need to pass stimulus function
         */        
        FitzHughNagumo1961OdeSystem fhn61_ode_system(pStimulus);
        
        /*
         * Create vectors of variable names & units
         * and set initial conditions
         */
         
        std::vector<std::string> variable_names;
        std::vector<std::string> variable_units;
        std::vector<double> initial_conditions;
        
        variable_names.push_back("V");
        variable_units.push_back("mV");
        initial_conditions.push_back(0.0);
        
        variable_names.push_back("w");
        variable_units.push_back("");
        initial_conditions.push_back(0.0);
        
        /*
         * Solve and write to file
         */
        runOdeSolverWithIonicModel(&fhn61_ode_system,
                                   endTime,
                                   timeStep,
                                   initial_conditions,
                                   pFilename,
                                   variable_names,
                                   variable_units); 
         
                                  
    }
    
    // Test Ode Solver for FHN61
    void testOdeSolverForFHN61WithInitialStimulus(void)
    {
        /*
         * Choose function for stimulus
         */             
        double magnitude = 1.0;  
        double duration  = 0.5;  // ms                   
        InitialStimulus stimulus(magnitude, duration); 
        
        /*
         * Solve 
         */
        runOdeSolverForFHN61(&stimulus, "FHN61Result.dat", 500.0, 0.01);
    }
    
    void testOdeSolverForFHN61WithRegularStimulus(void)
    {
        /*
         * Choose function for stimulus
         */             
        double magnitude = 1.0;  
        double duration  = 0.5 ;  // ms                     
        double frequency = 1.0/500.0;
        double start_stimulus = 40.0;                 
        RegularStimulus stimulus(magnitude, duration, frequency, start_stimulus); 

        /*
         * Solve 
         */
        runOdeSolverForFHN61(&stimulus, "FHN61RegResult.dat", 500.0, 0.01);
    }
    
        
    
    void runOdeSolverForLR91(AbstractStimulusFunction *pStimulus,
                             const char *pFilename,
                             double endTime,
                             double timeStep)
    {
        /*
         * Instantiate the ionic model: need to pass stimulus function
         */        
        LuoRudyIModel1991OdeSystem lr91_ode_system(pStimulus);
        
        /*
         * Create vectors of variable names & units
         */
        std::vector<std::string> variable_names;
        std::vector<std::string> variable_units;
        std::vector<double> initial_conditions;
        
        variable_names.push_back("h");
        variable_units.push_back("");
        initial_conditions.push_back(0.9833);
         
        variable_names.push_back("j");
        variable_units.push_back("");
        initial_conditions.push_back(0.9895);
         
        variable_names.push_back("m");
        variable_units.push_back("");
        initial_conditions.push_back(0.0017);
         
        variable_names.push_back("CaI");
        variable_units.push_back("mMol");
        initial_conditions.push_back(0.0002);
         
        variable_names.push_back("V");
        variable_units.push_back("mV");
        initial_conditions.push_back(-84.5);
         
        variable_names.push_back("d");
        variable_units.push_back("");
        initial_conditions.push_back(0.003);
         
        variable_names.push_back("f");
         variable_units.push_back("");
        initial_conditions.push_back(1);
         
        variable_names.push_back("x");
        variable_units.push_back("");
        initial_conditions.push_back(0.0056);
         
       
        /*
         * Solve and write to file
         */
        runOdeSolverWithIonicModel(&lr91_ode_system,
                                   endTime,
                                   timeStep,
                                   initial_conditions,
                                   pFilename,
                                   variable_names,
                                   variable_units);
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
        
        runOdeSolverForLR91(&stimulus, "NewDelayedStimLR91", end_time, time_step);                
        
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
                               
    }   
    
    
    
    
    

};



#endif //_TESTIONICMODELS_HPP_
