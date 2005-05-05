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
        double w = 0.0; // initial gating variable
        double magnitudeOfStimulus = -50.0;  
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
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus); 
        
        /*
         * Instantiate the LHN model: need to pass initial data and stimulus function
         */        
        AbstractOdeSystem *pFhn61OdeSystem = new FitzHughNagumo1961OdeSystem(pStimulus);
        
        /*
         * Choose an ode solver
         */      
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 500.0; // ms?
        double timeStep = 0.01;             
                
        OdeSolution solution = pMySolver->Solve(pFhn61OdeSystem, startTime, endTime, timeStep, initialConditions);  
        
              
        /*
         * Write data to a file FHN61.dat using ColumnDataWriter
         */                                                           
		        
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","FHN61");
        mpNewTestWriter->DefineFixedDimension("Time","ms", solution.mSolutions.size());
        int new_v_var_id = mpNewTestWriter->DefineVariable("V","mV");
        int new_time_var_id = mpNewTestWriter->DefineVariable("Time","ms");
        int new_w_var_id = mpNewTestWriter->DefineVariable("w"," ");
        mpNewTestWriter->EndDefineMode();
				
        for (int i = 0; i < solution.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(new_time_var_id, solution.mTime[i], i);
            mpNewTestWriter->PutVariable(new_v_var_id, solution.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(new_w_var_id, solution.mSolutions[i][1], i);            
        }
        mpNewTestWriter->Close();
        
        
//        //read in good data file and compare line by line
//        std::ifstream testfile("data/NewLR91.dat",std::ios::in);
//        std::ifstream goodfile("data/Lr91Good.dat",std::ios::in);
//        std::string teststring;
//        std::string goodstring;
//        while(getline(testfile, teststring))
//        {
//              getline(goodfile,goodstring);
//              TS_ASSERT_EQUALS(teststring,goodstring);
//        }
//        testfile.close();
//        goodfile.close();                              
    }
    
    // Test Ode Solver for FHN61
//    void testOdeSolverForFHN61WithRegularStimulus(void)
//    {
//        /*
//         * Set initial conditions and magnitude of stimulus
//         */
//		double voltage = 0.0; // initial resting potential
//        double w = 0.0; // initial gating variable
//        
//        double magnitudeOfStimulus = -50.0;  
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
//         * Write data to a file FHN61.dat using ColumnDataWriter
//         */                                                           
//		        
//        ColumnDataWriter *mpNewTestWriter;
//        mpNewTestWriter = new ColumnDataWriter("data","FHN61");
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
//        
//        
////        //read in good data file and compare line by line
////        std::ifstream testfile("data/NewLR91.dat",std::ios::in);
////        std::ifstream goodfile("data/Lr91Good.dat",std::ios::in);
////        std::string teststring;
////        std::string goodstring;
////        while(getline(testfile, teststring))
////        {
////              getline(goodfile,goodstring);
////              TS_ASSERT_EQUALS(teststring,goodstring);
////        }
////        testfile.close();
////        goodfile.close();                              
//    }	      
};

#endif //_TESTODESOLVERFORFHN61_HPP_
