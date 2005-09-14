#ifndef _TESTODESOLVERFORFHN61_HPP_
#define _TESTODESOLVERFORFHN61_HPP_

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
#include "FitzHughNagumo1961OdeSystem.hpp"


class TestOdeSolverForFHN61 : public CxxTest::TestSuite
{
    public:
    
    
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

#endif //_TESTODESOLVERFORFHN61_HPP_
