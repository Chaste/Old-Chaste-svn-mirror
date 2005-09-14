#ifndef _TESTODESOLVERFORLR91_HPP_
#define _TESTODESOLVERFORLR91_HPP_

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
#include "LuoRudyIModel1991OdeSystem.hpp"

/**
 * \todo Looking at the output files in Matlab, it seems to me that this model
 * is incorrect. Someone who knows more about it ought to take a look.
 *  - Jonathan C.
 */

class TestOdeSolverForLR91 : public CxxTest::TestSuite
{
    public:
    
    
    // Test Ode Solver for LR91
    void testOdeSolverForLR91WithInitialStimulus(void)
    {
        /*
         * Set initial conditions and magnitude of stimulus
         * 
         */
		double voltage = -84.5;
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;
        double d = 0.003;
        double f = 1;
        double x = 0.0056;
        double caI = 0.0002;
        double magnitudeOfStimulus = -80.0;  
        double durationOfStimulus  = 0.5 ;  // ms                     
        
        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initialConditions;
        initialConditions.push_back(h);
        initialConditions.push_back(j);
        initialConditions.push_back(m);
        initialConditions.push_back(caI);
        initialConditions.push_back(voltage);
        initialConditions.push_back(d);
        initialConditions.push_back(f);
        initialConditions.push_back(x);
         
        /*
         * Choose function for stimulus
         */             
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus); 
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        LuoRudyIModel1991OdeSystem lr91_ode_system(&stimulus);
        
        /*
         * Choose an ode solver
         */      
        EulerIvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 1.0;
        double timeStep = 0.01;             
                
        OdeSolution solution_new = solver.Solve(&lr91_ode_system, startTime, endTime, timeStep, initialConditions);
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
		        
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("testoutput","NewLR91");
        int new_time_var_id=mpNewTestWriter->DefineFixedDimension("Time","ms", solution_new.mSolutions.size());
        int new_v_var_id = mpNewTestWriter->DefineVariable("V","mV");
        int new_m_var_id = mpNewTestWriter->DefineVariable("m","");
        int new_h_var_id = mpNewTestWriter->DefineVariable("h","");
        int new_j_var_id = mpNewTestWriter->DefineVariable("j","");
        int new_d_var_id = mpNewTestWriter->DefineVariable("d","");
        int new_f_var_id = mpNewTestWriter->DefineVariable("f","");
        int new_x_var_id = mpNewTestWriter->DefineVariable("x","");
        int new_caI_var_id = mpNewTestWriter->DefineVariable("CaI","mMol");
        mpNewTestWriter->EndDefineMode();
				
        for (int i = 0; i < solution_new.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(new_time_var_id, solution_new.mTime[i], i);
            mpNewTestWriter->PutVariable(new_v_var_id, solution_new.mSolutions[i][4], i);
            mpNewTestWriter->PutVariable(new_m_var_id, solution_new.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(new_h_var_id, solution_new.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(new_j_var_id, solution_new.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(new_f_var_id, solution_new.mSolutions[i][6], i);
            mpNewTestWriter->PutVariable(new_d_var_id, solution_new.mSolutions[i][5], i);
            mpNewTestWriter->PutVariable(new_x_var_id, solution_new.mSolutions[i][7], i);
            mpNewTestWriter->PutVariable(new_caI_var_id, solution_new.mSolutions[i][3], i);            
        }
        
        delete mpNewTestWriter;
        
                               
    }	

	// Test Ode Solver for LR91
    void testOdeSolverForLR91NoStimulus(void)
    {
        /*
         * Set initial conditions and magnitude of stimulus
         * 
         */
        double voltage = -84.5;
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;
        double d = 0.003;
        double f = 1;
        double x = 0.0056;
        double caI = 0.0002;
        double magnitudeOfStimulus = 0.0;  
        double durationOfStimulus  = 0.5 ;  // ms                     
        
        /*
         * Collect initial data in a vector
         * 
         */   
        std::vector<double> initialConditions;
        initialConditions.push_back(h);
        initialConditions.push_back(j);
        initialConditions.push_back(m);
        initialConditions.push_back(caI);
        initialConditions.push_back(voltage);
        initialConditions.push_back(d);
        initialConditions.push_back(f);
        initialConditions.push_back(x);
         
        /*
         * Choose function for stimulus
         */             
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus); 
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        LuoRudyIModel1991OdeSystem lr91_ode_system(&stimulus);
        
        /*
         * Choose an ode solver
         */      
        EulerIvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 1.0;
        double timeStep = 0.01;             
                
        OdeSolution solution_new = solver.Solve(&lr91_ode_system, startTime, endTime, timeStep, initialConditions);
        
        
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
                
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("testoutput","NewNoStimLR91");
        int new_time_var_id=mpNewTestWriter->DefineFixedDimension("Time","ms", solution_new.mSolutions.size());
        int new_v_var_id = mpNewTestWriter->DefineVariable("V","mV");
        int new_m_var_id = mpNewTestWriter->DefineVariable("m","");
        int new_h_var_id = mpNewTestWriter->DefineVariable("h","");
        int new_j_var_id = mpNewTestWriter->DefineVariable("j","");
        int new_d_var_id = mpNewTestWriter->DefineVariable("d","");
        int new_f_var_id = mpNewTestWriter->DefineVariable("f","");
        int new_x_var_id = mpNewTestWriter->DefineVariable("x","");
        int new_caI_var_id = mpNewTestWriter->DefineVariable("CaI","mMol");
        mpNewTestWriter->EndDefineMode();
                
        for (int i = 0; i < solution_new.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(new_time_var_id, solution_new.mTime[i], i);
            mpNewTestWriter->PutVariable(new_v_var_id, solution_new.mSolutions[i][4], i);
            mpNewTestWriter->PutVariable(new_m_var_id, solution_new.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(new_h_var_id, solution_new.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(new_j_var_id, solution_new.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(new_f_var_id, solution_new.mSolutions[i][6], i);
            mpNewTestWriter->PutVariable(new_d_var_id, solution_new.mSolutions[i][5], i);
            mpNewTestWriter->PutVariable(new_x_var_id, solution_new.mSolutions[i][7], i);
            mpNewTestWriter->PutVariable(new_caI_var_id, solution_new.mSolutions[i][3], i);            
        }
        
        delete mpNewTestWriter;

                               
    }   
    
    // Test Ode Solver for LR91
    void testOdeSolverForLR91RegStimulus(void)
    {
        /*
         * This test doesn't really do anything. 
         * This is because it would have to run for quite a while and use up 
         * lots of memory to generate repeated stimulus results. To get lovely 
         * pictures increase the amount of time this runs for...
         * 
         */
         
         
        /*
         * Set initial conditions and magnitude of stimulus
         * 
         */
        double voltage = -84.5;
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;
        double d = 0.003;
        double f = 1;
        double x = 0.0056;
        double caI = 0.0002;

        double magnitudeOfStimulus = -80.0;  
        double durationOfStimulus  = 0.5 ;  // ms                     
        double frequency = 1.0/500.0;
        double startStimulus = 40.0;
        
        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initialConditions;
        initialConditions.push_back(h);
        initialConditions.push_back(j);
        initialConditions.push_back(m);
        initialConditions.push_back(caI);
        initialConditions.push_back(voltage);
        initialConditions.push_back(d);
        initialConditions.push_back(f);
        initialConditions.push_back(x);
        
        /*
         * Choose function for stimulus
         */             
        RegularStimulus stimulus(magnitudeOfStimulus, durationOfStimulus, frequency, startStimulus);
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        LuoRudyIModel1991OdeSystem lr91_ode_system(&stimulus);
        
        /*
         * Choose an ode solver
         */      
        EulerIvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 1.0;
        double timeStep = 0.01;             
                
        OdeSolution solution_new = solver.Solve(&lr91_ode_system, startTime, endTime, timeStep, initialConditions);  
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
                
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("testoutput","RegStimLR91");
        int new_time_var_id=mpNewTestWriter->DefineFixedDimension("Time","ms", solution_new.mSolutions.size());
        int new_v_var_id = mpNewTestWriter->DefineVariable("V","mV");
        int new_m_var_id = mpNewTestWriter->DefineVariable("m","");
        int new_h_var_id = mpNewTestWriter->DefineVariable("h","");
        int new_j_var_id = mpNewTestWriter->DefineVariable("j","");
        int new_d_var_id = mpNewTestWriter->DefineVariable("d","");
        int new_f_var_id = mpNewTestWriter->DefineVariable("f","");
        int new_x_var_id = mpNewTestWriter->DefineVariable("x","");
        int new_caI_var_id = mpNewTestWriter->DefineVariable("CaI","mMol");
        mpNewTestWriter->EndDefineMode();
                
        for (int i = 0; i < solution_new.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(new_time_var_id, solution_new.mTime[i], i);
            mpNewTestWriter->PutVariable(new_v_var_id, solution_new.mSolutions[i][4], i);
            mpNewTestWriter->PutVariable(new_m_var_id, solution_new.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(new_h_var_id, solution_new.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(new_j_var_id, solution_new.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(new_f_var_id, solution_new.mSolutions[i][6], i);
            mpNewTestWriter->PutVariable(new_d_var_id, solution_new.mSolutions[i][5], i);
            mpNewTestWriter->PutVariable(new_x_var_id, solution_new.mSolutions[i][7], i);
            mpNewTestWriter->PutVariable(new_caI_var_id, solution_new.mSolutions[i][3], i);            
        }
        
        delete mpNewTestWriter;
        
                                       
    } 
           
};

#endif //_TESTODESOLVERFORLR91_HPP_
