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
#include "TenTusscherModel2004OdeSystem.hpp"

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
   

    void runOdeSolverForHH52(AbstractStimulusFunction *pStimulus,
                             const char *pFilename,
                             double endTime,
                             double timeStep)
    {
        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(pStimulus);
        
        /*
         * Create vectors of variable names & units
         * and initial conditions
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
         * Set stimulus
         */
        double magnitude_of_stimulus = -20.0;  
        double duration_of_stimulus  = 0.5 ;  // ms                     
        InitialStimulus stimulus(magnitude_of_stimulus, duration_of_stimulus);

        runOdeSolverForHH52(&stimulus, "HH52Result.dat", 5.0, 0.0001);
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
         * Set stimulus
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
        LuoRudyIModel1991OdeSystem lr91_ode_system(pStimulus);
        
        /*
         * Create vectors of variable names & units
         * and set initial conditions
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




    void runOdeSolverForTT04(AbstractStimulusFunction *pStimulus,
                             const char *pFilename,
                             double endTime,
                             double timeStep)
    {
        TenTusscherModel2004OdeSystem tt04_ode_system(pStimulus);
        
        /*
         * Create vectors of variable names & units
         * and set initial conditions
         */
        std::vector<std::string> variable_names;
        std::vector<std::string> variable_units;
        std::vector<double> initial_conditions;

        variable_names.push_back("calcium_dynamics_Ca_i");
        variable_units.push_back("microMo");
        initial_conditions.push_back(0.0002);
        
        variable_names.push_back("calcium_dynamics_Ca_SR");
        variable_units.push_back("mMol");
        initial_conditions.push_back(0.2);
        
        variable_names.push_back("calcium_dynamics_g");
        variable_units.push_back("mMol");
        initial_conditions.push_back(1.0);
        
        variable_names.push_back("fast_sodium_current_h_gate_h");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(0.75);
        
        variable_names.push_back("fast_sodium_current_j_gate_j");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(0.75);
        
        variable_names.push_back("fast_sodium_current_m_gate_m");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(0.0);
        
        variable_names.push_back("L_type_calcium_current_d_gate_d");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(0.0);
        
        variable_names.push_back("L_type_calcium_current_f_Ca_gate_f_Ca");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(1.0);
        
        variable_names.push_back("L_type_calcium_current_f_gate_f");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(1.0);
        
        variable_names.push_back("membrane_V");
        variable_units.push_back("millivolts");
        initial_conditions.push_back(-86.2);
        
        variable_names.push_back("potassium_dynamics_K_i");
        variable_units.push_back("mMol");
        initial_conditions.push_back(138.3);
        
        variable_names.push_back("rapid_delayed_rectifier_current_X_r1_gate_X_r1");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(0.0);
        
        variable_names.push_back("rapid_delayed_rectifier_current_X_r2_gate_X_r2");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(1.0);
        
        variable_names.push_back("slow_delayed_rectifier_current_X_s_gate_X_s");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(0.0);
        
        variable_names.push_back("sodium_dynamics_Na_i");
        variable_units.push_back("mMol");
        initial_conditions.push_back(11.6);
        
        variable_names.push_back("transient_outward_current_r_gate_r");
        variable_units.push_back("milliamperes");
        initial_conditions.push_back(0.0);
        
        variable_names.push_back("transient_outward_current_s_gate_s");
        initial_conditions.push_back(1.0);
        variable_units.push_back("milliamperes");
        
        /*
         * Solve and write to file
         */
        runOdeSolverWithIonicModel(&tt04_ode_system,
                                   endTime,
                                   timeStep,
                                   initial_conditions,
                                   pFilename,
                                   variable_names,
                                   variable_units);
    }
    

    void TestTenTusscher04ModelWithNoStimulus(void)
    {
        /*
         * Set stimulus
         */
        double magnitude = 0.0;  
        double duration  = 1.0;  // ms                     
        InitialStimulus stimulus(magnitude, duration);
        
        double end_time = 10.0;  // ms
        double time_step = 0.01; // 1e-5 seconds in milliseconds           
        
        runOdeSolverForTT04(&stimulus, "NoStimulusTT04", end_time, time_step);  
    }

};



#endif //_TESTIONICMODELS_HPP_
