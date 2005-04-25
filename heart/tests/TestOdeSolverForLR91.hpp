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
#include "AbstractDataWriter.hpp"
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
        LR91Model *pLR91Model = new LR91Model(voltage, m, h, j, d, f, x, caI, pStimulus);
        
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
        
        OdeSolution solution = pLR91Model->SolveModel(startTime, endTime, timeStep,
                               intialConditions, pMySolver);
                            
        /*
         * Write data to a file LRresult.dat
         */          
        AbstractDataWriter *pResultWriter = new ColumnDataWriter(".", "LRResult");
        pResultWriter->DefineUnlimitedDimension("Time", "ms");
        int time_id = pResultWriter->DefineVariable("Time", "ms");
        int voltage_id = pResultWriter->DefineVariable("Voltage", "mV");
	    pResultWriter->EndDefineMode();
	    
	    for(int i=0;i<solution.mTime.size() ; i++)
	    {
	    	pResultWriter->PutVariable(time_id,solution.mTime[i]);
	    	pResultWriter->PutVariable(voltage_id,solution.mSolutions[i][0]);
			pResultWriter->AdvanceAlongUnlimitedDimension();
	    }
        pResultWriter->Close();                                                         
       // Solution.SaveToFile("LRresult.dat");
                               
    }

    // Test Ode Solver for LR91
    void testOdeSolverForLR91WithRegularStimulus(void)
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
        LR91Model *pLR91Model = new LR91Model(voltage, m, h, j, d, f, x, caI, pStimulus);
        
        /*
         * Choose an ode solver
         */      
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 1500.0;
        double timeStep = 0.001;             
        
        OdeSolution solution = pLR91Model->SolveModel(startTime, endTime, timeStep,
                               intialConditions, pMySolver);
                            
        /*
         * Write data to a file LRresult2.dat
         */                                  
                                  
         /*
         * Write data to a file LRresult.dat
         */          
//        AbstractDataWriter *pResultWriter = new ColumnDataWriter(".", "LRResult2");
//        pResultWriter->DefineUnlimitedDimension("Time", "ms");
//        int time_id = pResultWriter->DefineVariable("Time", "ms");
//        int voltage_id = pResultWriter->DefineVariable("Voltage", "mV");
//	    pResultWriter->EndDefineMode();
//	    
//	    for(int i=0;i<solution.mTime.size() ; i++)
//	    {
//	    	pResultWriter->PutVariable(time_id,solution.mTime[i]);
//	    	pResultWriter->PutVariable(voltage_id,solution.mSolutions[i][0]);
//			pResultWriter->AdvanceAlongUnlimitedDimension();
//	    }
//        pResultWriter->Close();             
                               
    }
//    // Test Ode Solver for LR91
//    void donttestOdeSolverForLR91WithRegularStimulus(void)
//    {
//        /*
//         * Set initial conditions and magnitude of stimulus
//         * 
//         */
//        double voltage = -75.0;
//        double m = 0.0017;
//        double h = 0.9833;
//        double j = 0.9895;
//        double d = 0.003;
//        double f = 1;
//        double x = 0.0056;
//        double caI = 0.0002;
//        double magnitudeOfStimulus = -80.0;  
//        double durationOfStimulus  = 0.5 ;  // ms
//        double frequencyOfStimulus = 1.0/300.0;
//        double startTimeOfStimulus = 50.0;
//        
//        /*
//         * Collect initiastl data in a vector
//         * 
//         */  
//        std::vector<double> intialConditions;
//        intialConditions.push_back(voltage);
//        intialConditions.push_back(m);
//        intialConditions.push_back(h);
//        intialConditions.push_back(j);
//        intialConditions.push_back(d);
//        intialConditions.push_back(f);
//        intialConditions.push_back(x);
//        intialConditions.push_back(caI);
//                     
//        /*
//         * Choose function for stimulus
//         */             
//        AbstractStimulusFunction *pStimulus = new RegularStimulus(magnitudeOfStimulus, durationOfStimulus,
//                                                                 frequencyOfStimulus, startTimeOfStimulus); 
//        
//        /*
//         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
//         */
//        LR91Model *pLR91Model = new LR91Model(voltage, m, h, j, d, f, x, caI, pStimulus);
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
//        double endTime = 600.0;
//        double timeStep = 0.01;             
//        
//        OdeSolution Solution = pLR91Model->SolveModel(startTime, endTime, timeStep,
//                               intialConditions, pMySolver);
//                            
//        /*
//         * Write data to a file LRresult2.dat
//         */                                  
//                                  
//      //  Solution.SaveToFile("LRresult2.dat");
//                               
//    }
	
     // Test Ode Solver for LR91
    void testOdeSolverForLR91WithRegularStimulusAndRungeKutta(void)
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
        LR91Model *pLR91Model = new LR91Model(voltage, m, h, j, d, f, x, caI, pStimulus);
        
        /*
         * Choose an ode solver
         */      
        AbstractIvpOdeSolver *pMySolver = new RungeKutta4IvpOdeSolver();
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 1500.0;
        double timeStep = 0.01;             
        
//        OdeSolution solution = pLR91Model->SolveModel(startTime, endTime, timeStep,
//                               intialConditions, pMySolver);
                            
        /*
         * Write data to a file LRresult3.dat
         */                                         
//        AbstractDataWriter *pResultWriter = new ColumnDataWriter(".", "LRResult3");
//        pResultWriter->DefineUnlimitedDimension("Time", "ms");
//        int time_id = pResultWriter->DefineVariable("Time", "ms");
//        int voltage_id = pResultWriter->DefineVariable("Voltage", "mV");
//	    pResultWriter->EndDefineMode();
//	    
//	    for(int i=0;i<solution.mTime.size() ; i++)
//	    {
//	    	pResultWriter->PutVariable(time_id,solution.mTime[i]);
//	    	pResultWriter->PutVariable(voltage_id,solution.mSolutions[i][0]);
//			pResultWriter->AdvanceAlongUnlimitedDimension();
//	    }
//        pResultWriter->Close();             
                               
    }
           
};



#endif //_TESTODESOLVERFORLR91_HPP_
