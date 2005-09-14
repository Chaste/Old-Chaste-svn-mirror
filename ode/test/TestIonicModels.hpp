#ifndef _TESTIONICMODELS_HPP_
#define _TESTIONICMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

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

#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"

const double TOLERANCE = 1e-2; // Not used at present


class TestOdeSolverForHH52 : public CxxTest::TestSuite
{
    public:
    
    
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
         * Set initial conditions
         * 
         */
         
        double voltage = 0.0;
        double n = 0.31768;
        double h = 0.59612;
        double m = 0.05293;

        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initial_conditions;
        initial_conditions.push_back(voltage);
        initial_conditions.push_back(n);
        initial_conditions.push_back(h);
        initial_conditions.push_back(m);        

        /*
         * Instantiate the ionic model: need to pass stimulus function
         */        
        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(pStimulus);
        
        /*
         * Choose an ode solver
         */
        RungeKutta4IvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double start_time = 0.0;
        int step_per_row = 10;             
                
        OdeSolution solution_new = solver.Solve(&hh52_ode_system, start_time, endTime, timeStep, initial_conditions);
        
              
        /*
         * Write data to a file using ColumnDataWriter
         */                                                           
        ColumnDataWriter *p_writer;
        p_writer = new ColumnDataWriter("testoutput",pFilename);
        int time_var_id=p_writer->DefineUnlimitedDimension("Time","ms"); //, solution_new.mSolutions.size());
        int v_var_id = p_writer->DefineVariable("V","milliamperes");
        int n_var_id = p_writer->DefineVariable("n","");
        int h_var_id = p_writer->DefineVariable("h","");
        int m_var_id = p_writer->DefineVariable("m","");
        p_writer->EndDefineMode();
                
        for (int i = 0; i < solution_new.mSolutions.size(); i+=step_per_row) 
        {
            if (i!=0)
            {
                p_writer->AdvanceAlongUnlimitedDimension();
            }
            p_writer->PutVariable(time_var_id, solution_new.mTime[i]);
            p_writer->PutVariable(v_var_id, solution_new.mSolutions[i][0]);
            p_writer->PutVariable(n_var_id, solution_new.mSolutions[i][1]);
            p_writer->PutVariable(h_var_id, solution_new.mSolutions[i][2]);
            p_writer->PutVariable(m_var_id, solution_new.mSolutions[i][3]);
            
        }
        
        delete p_writer;

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
         * This test doesn't really do anything. 
         * This is because it would have to run for quite a while and use up 
         * lots of memory to generate repeated stimulus results. To get lovely 
         * pictures increase the amount of time this runs for...
         * 
         */
        
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
         * Set initial conditions
         */
        double voltage = 0.0; // initial resting potential
        double w = 0.0;       // initial value for recovery variable
        
        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initial_conditions;
        initial_conditions.push_back(voltage);
        initial_conditions.push_back(w);
        
        /*
         * Instantiate the FHN model: need to pass stimulus function
         */        
        FitzHughNagumo1961OdeSystem fhn61_ode_system(pStimulus);
        
        /*
         * Choose an ode solver
         */      
        EulerIvpOdeSolver solver;
        
        /*
         * Solve
         */
        double start_time = 0.0;
        OdeSolution solution = solver.Solve(&fhn61_ode_system,
                                            start_time,
                                            endTime,
                                            timeStep,
                                            initial_conditions);  
         
        /*
         * Write data to a file FHN61.dat using ColumnDataWriter
         */                                                           
        ColumnDataWriter writer("testoutput", pFilename);
        int new_time_var_id = writer.DefineFixedDimension("Time", "ms",
                                                          solution.mSolutions.size());
        int new_v_var_id = writer.DefineVariable("V","mV");
        int new_w_var_id = writer.DefineVariable("w","");
        writer.EndDefineMode();
                
        for (int i = 0; i < solution.mSolutions.size(); i++) 
        {
            writer.PutVariable(new_time_var_id, solution.mTime[i], i);
            writer.PutVariable(new_v_var_id, solution.mSolutions[i][0], i);
            writer.PutVariable(new_w_var_id, solution.mSolutions[i][1], i);            
        }
        writer.Close();                            
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
};



#endif //_TESTIONICMODELS_HPP_
