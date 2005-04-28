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
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus); 
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        AbstractOdeSystem *pLr91OdeSystem = new LuoRudyIModel1991OdeSystem(pStimulus);
        
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
                
        OdeSolution SolutionNew = pMySolver->Solve(pLr91OdeSystem, startTime, endTime, timeStep, initialConditions);  
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
		        
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","NewLR91");
        mpNewTestWriter->DefineFixedDimension("Time","ms", SolutionNew.mSolutions.size());
        int new_v_var_id = mpNewTestWriter->DefineVariable("V","milliamperes");
        int new_time_var_id = mpNewTestWriter->DefineVariable("Time","ms");
        int new_m_var_id = mpNewTestWriter->DefineVariable("m"," ");
        int new_h_var_id = mpNewTestWriter->DefineVariable("h"," ");
        int new_j_var_id = mpNewTestWriter->DefineVariable("j"," ");
        int new_d_var_id = mpNewTestWriter->DefineVariable("d"," ");
        int new_f_var_id = mpNewTestWriter->DefineVariable("f"," ");
        int new_x_var_id = mpNewTestWriter->DefineVariable("x"," ");
        int new_caI_var_id = mpNewTestWriter->DefineVariable("CaI","mMol ");
        mpNewTestWriter->EndDefineMode();
				
        for (int i = 0; i < SolutionNew.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(new_time_var_id, SolutionNew.mTime[i], i);
            mpNewTestWriter->PutVariable(new_v_var_id, SolutionNew.mSolutions[i][4], i);
            mpNewTestWriter->PutVariable(new_m_var_id, SolutionNew.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(new_h_var_id, SolutionNew.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(new_j_var_id, SolutionNew.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(new_f_var_id, SolutionNew.mSolutions[i][6], i);
            mpNewTestWriter->PutVariable(new_d_var_id, SolutionNew.mSolutions[i][5], i);
            mpNewTestWriter->PutVariable(new_x_var_id, SolutionNew.mSolutions[i][7], i);
            mpNewTestWriter->PutVariable(new_caI_var_id, SolutionNew.mSolutions[i][3], i);            
        }
        mpNewTestWriter->Close();
        
        
        //read in good data file and compare line by line
        std::ifstream testfile("data/NewLR91.dat",std::ios::in);
        std::ifstream goodfile("data/Lr91Good.dat",std::ios::in);
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
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus); 
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        AbstractOdeSystem *pLr91OdeSystem = new LuoRudyIModel1991OdeSystem(pStimulus);
        
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
                
        OdeSolution SolutionNew = pMySolver->Solve(pLr91OdeSystem, startTime, endTime, timeStep, initialConditions);  
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
                
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","NewNoStimLR91");
        mpNewTestWriter->DefineFixedDimension("Time","ms", SolutionNew.mSolutions.size());
        int new_v_var_id = mpNewTestWriter->DefineVariable("V","milliamperes");
        int new_time_var_id = mpNewTestWriter->DefineVariable("Time","ms");
        int new_m_var_id = mpNewTestWriter->DefineVariable("m"," ");
        int new_h_var_id = mpNewTestWriter->DefineVariable("h"," ");
        int new_j_var_id = mpNewTestWriter->DefineVariable("j"," ");
        int new_d_var_id = mpNewTestWriter->DefineVariable("d"," ");
        int new_f_var_id = mpNewTestWriter->DefineVariable("f"," ");
        int new_x_var_id = mpNewTestWriter->DefineVariable("x"," ");
        int new_caI_var_id = mpNewTestWriter->DefineVariable("CaI","mMol ");
        mpNewTestWriter->EndDefineMode();
                
        for (int i = 0; i < SolutionNew.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(new_time_var_id, SolutionNew.mTime[i], i);
            mpNewTestWriter->PutVariable(new_v_var_id, SolutionNew.mSolutions[i][4], i);
            mpNewTestWriter->PutVariable(new_m_var_id, SolutionNew.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(new_h_var_id, SolutionNew.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(new_j_var_id, SolutionNew.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(new_f_var_id, SolutionNew.mSolutions[i][6], i);
            mpNewTestWriter->PutVariable(new_d_var_id, SolutionNew.mSolutions[i][5], i);
            mpNewTestWriter->PutVariable(new_x_var_id, SolutionNew.mSolutions[i][7], i);
            mpNewTestWriter->PutVariable(new_caI_var_id, SolutionNew.mSolutions[i][3], i);            
        }
        mpNewTestWriter->Close();
        
        
        //read in good data file and compare line by line
        std::ifstream testfile("data/NewNoStimLR91.dat",std::ios::in);
        std::ifstream goodfile("data/Lr91NoStimGood.dat",std::ios::in);
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
        AbstractStimulusFunction *pStimulus = new RegularStimulus(magnitudeOfStimulus, durationOfStimulus, frequency, startStimulus); 
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        AbstractOdeSystem *pLr91OdeSystem = new LuoRudyIModel1991OdeSystem(pStimulus);
        
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
                
        OdeSolution SolutionNew = pMySolver->Solve(pLr91OdeSystem, startTime, endTime, timeStep, initialConditions);  
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
                
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","RegStimLR91");
        mpNewTestWriter->DefineFixedDimension("Time","ms", SolutionNew.mSolutions.size());
        int new_v_var_id = mpNewTestWriter->DefineVariable("V","milliamperes");
        int new_time_var_id = mpNewTestWriter->DefineVariable("Time","ms");
        int new_m_var_id = mpNewTestWriter->DefineVariable("m"," ");
        int new_h_var_id = mpNewTestWriter->DefineVariable("h"," ");
        int new_j_var_id = mpNewTestWriter->DefineVariable("j"," ");
        int new_d_var_id = mpNewTestWriter->DefineVariable("d"," ");
        int new_f_var_id = mpNewTestWriter->DefineVariable("f"," ");
        int new_x_var_id = mpNewTestWriter->DefineVariable("x"," ");
        int new_caI_var_id = mpNewTestWriter->DefineVariable("CaI","mMol ");
        mpNewTestWriter->EndDefineMode();
                
        for (int i = 0; i < SolutionNew.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(new_time_var_id, SolutionNew.mTime[i], i);
            mpNewTestWriter->PutVariable(new_v_var_id, SolutionNew.mSolutions[i][4], i);
            mpNewTestWriter->PutVariable(new_m_var_id, SolutionNew.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(new_h_var_id, SolutionNew.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(new_j_var_id, SolutionNew.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(new_f_var_id, SolutionNew.mSolutions[i][6], i);
            mpNewTestWriter->PutVariable(new_d_var_id, SolutionNew.mSolutions[i][5], i);
            mpNewTestWriter->PutVariable(new_x_var_id, SolutionNew.mSolutions[i][7], i);
            mpNewTestWriter->PutVariable(new_caI_var_id, SolutionNew.mSolutions[i][3], i);            
        }
        mpNewTestWriter->Close();
        
        
//        //read in good data file and compare line by line
//        std::ifstream testfile("data/NewNoStimLR91.dat",std::ios::in);
//        std::ifstream goodfile("data/Lr91NoStimGood.dat",std::ios::in);
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



#endif //_TESTODESOLVERFORLR91_HPP_
