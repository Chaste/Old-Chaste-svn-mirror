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
    
    
    // Test Ode Solver for FHN61
    void testOdeSolverForFHN61WithInitialStimulus(void)
    {
        /*
         * Set initial conditions and magnitude of stimulus
         */
		double voltage = 0.0; // initial resting potential
        double w = 0.0; // initial value for gating variable
        double magnitudeOfStimulus = 1.0;  
        double durationOfStimulus  = 0.5;  // ms                   
        
        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initialConditions;
        initialConditions.push_back(voltage);
        initialConditions.push_back(w);
         
        /*
         * Choose function for stimulus
         */             
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus); 
        
        /*
         * Instantiate the LHN model: need to pass initial data and stimulus function
         */        
        FitzHughNagumo1961OdeSystem Fhn61OdeSystem(&stimulus);
        
        /*
         * Choose an ode solver
         */      
        EulerIvpOdeSolver mySolver;
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 500.0; // ms
        double timeStep = 0.01;             
                
        OdeSolution solution = mySolver.Solve(&Fhn61OdeSystem, startTime, endTime, timeStep, initialConditions);  
        
              
        /*
         * Write data to a file FHN61.dat using ColumnDataWriter
         */                                                           
		        
        ColumnDataWriter mNewTestWriter("testoutput","FHN61");
        int new_time_var_id=mNewTestWriter.DefineFixedDimension("Time","ms", solution.mSolutions.size());
        int new_v_var_id = mNewTestWriter.DefineVariable("V","mV");
        int new_w_var_id = mNewTestWriter.DefineVariable("w"," ");
        mNewTestWriter.EndDefineMode();
				
        for (int i = 0; i < solution.mSolutions.size(); i++) 
        {
            mNewTestWriter.PutVariable(new_time_var_id, solution.mTime[i], i);
            mNewTestWriter.PutVariable(new_v_var_id, solution.mSolutions[i][0], i);
            mNewTestWriter.PutVariable(new_w_var_id, solution.mSolutions[i][1], i);            
        }
        mNewTestWriter.Close();                            
    }
    
//    void testOdeSolverForFHN61WithRegularStimulus(void)
//    {
//        /*
//         * Set initial conditions and magnitude of stimulus
//         */
//		double voltage = 0.0; // initial resting potential
//        double w = 0.0; // initial gating variable
//        
//        double magnitudeOfStimulus = 1.0;  
//        double durationOfStimulus  = 0.5 ;  // ms                     
//        double frequency = 1.0/500.0;
//        double startStimulus = 40.0;                 
//        
//        /*
//         * Collect initial data in a vector
//         * 
//         */  
//        std::vector<double> initialConditions;
//        initialConditions.push_back(voltage);
//        initialConditions.push_back(w);
//         
//        /*
//         * Choose function for stimulus
//         */             
//        AbstractStimulusFunction *pStimulus = new RegularStimulus(magnitudeOfStimulus, durationOfStimulus, frequency, startStimulus); 
//        
//        /*
//         * Instantiate the LHN model: need to pass initial data and stimulus function
//         */        
//        AbstractOdeSystem *pFhn61OdeSystem = new FitzHughNagumo1961OdeSystem(pStimulus);
//        
//        /*
//         * Choose an ode solver
//         */      
//        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
//        
//        /*
//         * Solve 
//         */
//        double startTime = 0.0;
//        double endTime = 500.0; // ms?
//        double timeStep = 0.01;             
//                
//        OdeSolution solution = pMySolver->Solve(pFhn61OdeSystem, startTime, endTime, timeStep, initialConditions);  
//        
//              
//        /*
//         * Write data to a file FHN61a.dat using ColumnDataWriter
//         */                                                           
//		        
//        ColumnDataWriter *mpNewTestWriter;
//        mpNewTestWriter = new ColumnDataWriter("testoutput","FHN61a");
//        mpNewTestWriter->DefineFixedDimension("Time","ms", solution.mSolutions.size());
//        int new_v_var_id = mpNewTestWriter->DefineVariable("V","mV");
//        int new_time_var_id = mpNewTestWriter->DefineVariable("Time","ms");
//        int new_w_var_id = mpNewTestWriter->DefineVariable("w"," ");
//        mpNewTestWriter->EndDefineMode();
//				
//        for (int i = 0; i < solution.mSolutions.size(); i++) 
//        {
//            mpNewTestWriter->PutVariable(new_time_var_id, solution.mTime[i], i);
//            mpNewTestWriter->PutVariable(new_v_var_id, solution.mSolutions[i][0], i);
//            mpNewTestWriter->PutVariable(new_w_var_id, solution.mSolutions[i][1], i);            
//        }
//        mpNewTestWriter->Close();                             
//    }	      
};

#endif //_TESTODESOLVERFORFHN61_HPP_
