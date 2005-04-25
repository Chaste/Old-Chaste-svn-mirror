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
#include "SodiumCurrentLR91.hpp"
#include "PotassiumTimeDependentCurrentLR91.hpp"
#include "PotassiumTimeIndependentCurrentLR91.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"


class TestOdeSolverForLR91 : public CxxTest::TestSuite
{
    public:
    
    // Test Ode Solver for LR91
    void testOdeSolverForLR91(void)
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
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus); 
        
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
        
        OdeSolution Solution = pLR91Model->SolveModel(startTime, endTime, timeStep,
                               intialConditions, pMySolver);
                            
        /*
         * Write data to a file LRresult.dat
         */                                                           
        Solution.SaveToFile("LRresult.dat");
                               
    }

	// Tests that Na gating variables within range 0 to 1 inclusive
	void testSodiumGatingVariables( void )
	{
		double m = 0.0017;
		double h = 0.9833;
		double j = 0.9895;
		double voltage = -75.0;
		
		SodiumCurrentLR91 *pINa = new SodiumCurrentLR91();
		pINa->UpdateMagnitudeOfCurrent(voltage,m,h,j);
		
		std::cout << "\n";
		std::cout << "INa m gate: " << pINa->GetM() << "\n";
		TS_ASSERT(pINa->GetM() >= 0);
		TS_ASSERT(pINa->GetM() <= 1);
		
		std::cout << "INa h gate: " << pINa->GetH() << "\n";
		TS_ASSERT(pINa->GetH() >= 0);
		TS_ASSERT(pINa->GetH() <= 1);		
	}
	
	// Tests that K gating variables within range 0 to 1 inclusive
	void testPotassiumTimeDependentGatingVariables( void )
	{
		double x = 0.0056;
		double voltage = -75.0;
		
		PotassiumTimeDependentCurrentLR91 *pIK = new PotassiumTimeDependentCurrentLR91();
		pIK->UpdateMagnitudeOfCurrent(voltage,x);
			
		std::cout << "IK n gate: " << pIK->GetX() << "\n";
		TS_ASSERT(pIK->GetX() >= 0);
		TS_ASSERT(pIK->GetX() <= 1);		
	}
	
	// Tests that K gating variables within range 0 to 1 inclusive
	void testPotassiumTimeIndependentGatingVariables( void )
	{
		double voltage = -75.0;
		
		PotassiumTimeIndependentCurrentLR91 *pIK = new PotassiumTimeIndependentCurrentLR91();
		pIK->UpdateMagnitudeOfCurrent(voltage);
			
		std::cout << "IK n gate: " << pIK->GetK1(voltage) << "\n";
		TS_ASSERT(pIK->GetK1(voltage) >= 0);
		TS_ASSERT(pIK->GetK1(voltage) <= 1);		
	}
	
	// Tests that K gating variables within range 0 to 1 inclusive
	void testPlateauPotassiumCurrentLR91( void )
	{
		double x = 0.0056;
		double voltage = -75.0;
		
		PotassiumTimeDependentCurrentLR91 *pIK = new PotassiumTimeDependentCurrentLR91();
		pIK->UpdateMagnitudeOfCurrent(voltage,x);
			
		std::cout << "IK n gate: " << pIK->GetX() << "\n";
		TS_ASSERT(pIK->GetX() >= 0);
		TS_ASSERT(pIK->GetX() <= 1);		
	}
	
	
	
    
};



#endif //_TESTODESOLVERFORLR91_HPP_
