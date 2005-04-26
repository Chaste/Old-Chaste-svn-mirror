#ifndef _TESTODESOLVERFORLR91_HPP_
#define _TESTODESOLVERFORLR91_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "LR91Model.hpp"
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"
#include "SodiumCurrentLR91.hpp"
#include "PotassiumTimeDependentCurrentLR91.hpp"
#include "PotassiumTimeIndependentCurrentLR91.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AdamsBashforthIvpOdeSolver.hpp"
#include "OdeSolution.hpp"
#include "ColumnDataWriter.hpp"


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
		double voltage = -75.0;
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
        std::vector<double> intialConditions;
        intialConditions.push_back(voltage);
        intialConditions.push_back(m);
        intialConditions.push_back(h);
        intialConditions.push_back(j);
        intialConditions.push_back(d);
        intialConditions.push_back(f);
        intialConditions.push_back(x);
        intialConditions.push_back(caI);
		             
        /*
         * Choose function for stimulus
         */             
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus); 
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */
        LR91Model *pLR91Model = new LR91Model(pStimulus);
        
        /*
         * Choose an ode solver
         */      
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 600.0;
        double timeStep = 0.01;             
        
        OdeSolution Solution = pLR91Model->SolveModel(startTime, endTime, timeStep,
                               intialConditions, pMySolver);
        
        // std::cout << "solved the system" << std::endl;
                               
        // Test that voltage is within [eK, eNa] range, given Istim = -80
        for (int i = 0; i < Solution.mSolutions.size(); i++)                   
        {
            TS_ASSERT(eK < Solution.mSolutions[i][0]);
            TS_ASSERT(eNa > Solution.mSolutions[i][0]);
        }
     
     
//        /*
//         * Write data to a file testLR91.dat using ColumnDataWriter
//         */                                                           
//		
//		ColumnDataWriter *mpTestWriter;
//		mpTestWriter = new ColumnDataWriter("data","testLR91");
//		mpTestWriter->DefineFixedDimension("Time","ms", Solution.mSolutions.size());
//		int v_var_id = mpTestWriter->DefineVariable("V","milliamperes");
//		int m_var_id = mpTestWriter->DefineVariable("m"," ");
//		int h_var_id = mpTestWriter->DefineVariable("h"," ");
//		int j_var_id = mpTestWriter->DefineVariable("j"," ");
//		int d_var_id = mpTestWriter->DefineVariable("d"," ");
//		int f_var_id = mpTestWriter->DefineVariable("f"," ");
//		int x_var_id = mpTestWriter->DefineVariable("x"," ");
//		int caI_var_id = mpTestWriter->DefineVariable("CaI","mMol ");
//		mpTestWriter->EndDefineMode();
//		
//		
//		for (int i = 0; i < Solution.mSolutions.size(); i++) 
//		{
//		  	mpTestWriter->PutVariable(v_var_id, Solution.mSolutions[i][0], i);
//		  	mpTestWriter->PutVariable(m_var_id, Solution.mSolutions[i][1], i);
//		  	mpTestWriter->PutVariable(h_var_id, Solution.mSolutions[i][2], i);
//		  	mpTestWriter->PutVariable(j_var_id, Solution.mSolutions[i][3], i);
//		  	mpTestWriter->PutVariable(f_var_id, Solution.mSolutions[i][4], i);
//		  	mpTestWriter->PutVariable(d_var_id, Solution.mSolutions[i][5], i);
//		  	mpTestWriter->PutVariable(x_var_id, Solution.mSolutions[i][6], i);
//		  	mpTestWriter->PutVariable(caI_var_id, Solution.mSolutions[i][7], i);        	
//		}
//		mpTestWriter->Close();
		
                               
    }
	
     // Test Ode Solver for LR91
    void donottestOdeSolverForLR91WithRegularStimulusAndRungeKutta(void)
    {
        /*
         * Set initial conditions and magnitude of stimulus
         * 
         */
        double voltage = -75.0;
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;
        double d = 0.003;
        double f = 1;
        double x = 0.0056;
        double caI = 0.0002;
        double magnitudeOfStimulus = -80.0;  
        double durationOfStimulus  = 0.5 ;  // ms
        double frequencyOfStimulus = 1.0/300.0;
        double startTimeOfStimulus = 50.0;
        
        /*
         * Collect initiastl data in a vector
         * 
         */  
        std::vector<double> intialConditions;
        intialConditions.push_back(voltage);
        intialConditions.push_back(m);
        intialConditions.push_back(h);
        intialConditions.push_back(j);
        intialConditions.push_back(d);
        intialConditions.push_back(f);
        intialConditions.push_back(x);
        intialConditions.push_back(caI);
                     
        /*
         * Choose function for stimulus
         */             
        AbstractStimulusFunction *pStimulus = new RegularStimulus(magnitudeOfStimulus, durationOfStimulus,
                                                                 frequencyOfStimulus, startTimeOfStimulus); 
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */
        LR91Model *pLR91Model = new LR91Model(pStimulus);
        
        /*
         * Choose an ode solver
         */      
        AbstractIvpOdeSolver *pMySolver = new RungeKutta4IvpOdeSolver();
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 600.0;
        double timeStep = 0.01;             
        
        OdeSolution Solution = pLR91Model->SolveModel(startTime, endTime, timeStep,
                               intialConditions, pMySolver);
                            
        /*
         * Write data to a file LRresult3.dat
         */                                  
                                  
        ColumnDataWriter *myDataWriter=new ColumnDataWriter("data","LRresult3");
		myDataWriter->DefineFixedDimension("Time","mS",Solution.mSolutions.size());
		int timeVariable = myDataWriter->DefineVariable("Time","msecs");
		int voltageVariable = myDataWriter->DefineVariable("Voltage","V");
		myDataWriter->EndDefineMode();
		for (int i=0; i<Solution.mSolutions.size(); i++)
		{
			myDataWriter->PutVariable(timeVariable,Solution.mTime[i],i);	
			myDataWriter->PutVariable(voltageVariable,Solution.mSolutions[i][0],i);
		}
		myDataWriter->Close();                           
    }
           
};



#endif //_TESTODESOLVERFORLR91_HPP_
