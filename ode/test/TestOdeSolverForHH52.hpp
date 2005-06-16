#ifndef _TESTODESOLVERFORHH52_HPP_
#define _TESTODESOLVERFORHH52_HPP_

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
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"


class TestOdeSolverForHH52 : public CxxTest::TestSuite
{
    public:
    
    
    // Test Ode Solver for HH52
    void testOdeSolverForHH52WithInitialStimulus(void)
    {
        /*
         * Set initial conditions and magnitude of stimulus
         * 
         */
		double voltage = 0.0;
        double n = 0.31768;
        double h = 0.59612;
        double m = 0.05293;
  
        double magnitudeOfStimulus = -20.0;  
        double durationOfStimulus  = 0.5 ;  // ms                     
        
        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initialConditions;
        initialConditions.push_back(voltage);
        initialConditions.push_back(n);
        initialConditions.push_back(h);
        initialConditions.push_back(m);        
        /*
         * Choose function for stimulus
         */             
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus); 
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        AbstractOdeSystem *pHh52OdeSystem = new HodgkinHuxleySquidAxon1952OriginalOdeSystem(pStimulus);
        
        /*
         * Choose an ode solver
         */      
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 1.0;
        double timeStep = 0.01;             
                
        OdeSolution SolutionNew = pMySolver->Solve(pHh52OdeSystem, startTime, endTime, timeStep, initialConditions);  
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
		        
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","HH52Result");
        mpNewTestWriter->DefineFixedDimension("Time","ms", SolutionNew.mSolutions.size());
        int time_var_id = mpNewTestWriter->DefineVariable("Time","ms");
        int v_var_id = mpNewTestWriter->DefineVariable("V","milliamperes");
        int n_var_id = mpNewTestWriter->DefineVariable("n"," ");
        int h_var_id = mpNewTestWriter->DefineVariable("h"," ");
        int m_var_id = mpNewTestWriter->DefineVariable("m"," ");
        mpNewTestWriter->EndDefineMode();
				
        for (int i = 0; i < SolutionNew.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(time_var_id, SolutionNew.mTime[i], i);
            mpNewTestWriter->PutVariable(v_var_id, SolutionNew.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(n_var_id, SolutionNew.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(h_var_id, SolutionNew.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(m_var_id, SolutionNew.mSolutions[i][3], i);         
        }
        mpNewTestWriter->Close();
        
        
        //read in good data file and compare line by line
        std::ifstream testfile("data/HH52Result.dat",std::ios::in);
        std::ifstream goodfile("data/HH52ResultGood.dat",std::ios::in);
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
         * Set initial conditions and magnitude of stimulus
         * 
         */
        double voltage = 0.0;
        double n = 0.31768;
        double h = 0.59612;
        double m = 0.05293;
  
        double magnitudeOfStimulus = -20.0;  
        double durationOfStimulus  = 0.5 ;  // ms                     
        double frequency = 1.0/50.0;
        double startStimulus = 40.0;              
        
        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initialConditions;
        initialConditions.push_back(voltage);
        initialConditions.push_back(n);
        initialConditions.push_back(h);
        initialConditions.push_back(m);        
        /*
         * Choose function for stimulus
         */             
       
        AbstractStimulusFunction *pStimulus = new RegularStimulus(magnitudeOfStimulus, durationOfStimulus, frequency, startStimulus); 
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        AbstractOdeSystem *pHh52OdeSystem = new HodgkinHuxleySquidAxon1952OriginalOdeSystem(pStimulus);
        
        /*
         * Choose an ode solver
         */      
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 1.0;
        double timeStep = 0.01;             
                
        OdeSolution SolutionNew = pMySolver->Solve(pHh52OdeSystem, startTime, endTime, timeStep, initialConditions);  
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
                
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","HH52RegResult");
        mpNewTestWriter->DefineFixedDimension("Time","ms", SolutionNew.mSolutions.size());
        int time_var_id = mpNewTestWriter->DefineVariable("Time","ms");
        int v_var_id = mpNewTestWriter->DefineVariable("V","milliamperes");
        int n_var_id = mpNewTestWriter->DefineVariable("n"," ");
        int h_var_id = mpNewTestWriter->DefineVariable("h"," ");
        int m_var_id = mpNewTestWriter->DefineVariable("m"," ");
        mpNewTestWriter->EndDefineMode();
                
        for (int i = 0; i < SolutionNew.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(time_var_id, SolutionNew.mTime[i], i);
            mpNewTestWriter->PutVariable(v_var_id, SolutionNew.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(n_var_id, SolutionNew.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(h_var_id, SolutionNew.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(m_var_id, SolutionNew.mSolutions[i][3], i);         
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

};



#endif //_TESTODESOLVERFORHH52_HPP_
