#ifndef _TESTLR91CURRENTS_HPP_
#define _TESTLR91CURRENTS_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include "ConstantsLR91.hpp"
#include "SodiumCurrentLR91.hpp"
#include "SlowInwardCurrentLR91.hpp"
#include "PotassiumTimeIndependentCurrentLR91.hpp"
#include "PotassiumTimeDependentCurrentLR91.hpp"
#include "PlateauPotassiumCurrentLR91.hpp"
#include "BackgroundCurrentLR91.hpp"

class TestLR91Currents : public CxxTest::TestSuite
{
	public:
	
	// Test fast sodium current, INa
	void TestSodiumCurrentLR91(void)
	{
		SodiumCurrentLR91  *mySodiumCurrent ;
		mySodiumCurrent =  new SodiumCurrentLR91 () ;
	}
	
	// Test slow inward current, ISi
	void TestSlowInwardCurrentLR91(void)
	{
		SlowInwardCurrentLR91  *mySlowInwardCurrent;
		mySlowInwardCurrent =  new SlowInwardCurrentLR91();
	}
//	
//	// Test potassium time dependent current, IK
//	void dontTestPotassiumTimeDependentCurrentLR91(void)
//	{   
//		PotassiumTimeDependentCurrentLR91 *myPotassiumTimeDependentCurrent;
//		myPotassiumTimeDependentCurrent = new PotassiumTimeDependentCurrentLR91();
//		double v = 20.0;
//		myPotassiumTimeDependentCurrent->SetGatingVariables(v);
//		double xtest=myPotassiumTimeDependentCurrent->GetX();
//		double xitest =myPotassiumTimeDependentCurrent->GetXi(v);
//		TS_ASSERT_DELTA(myPotassiumTimeDependentCurrent->GetMagnitudeOfCurrent() , gK*xtest*xitest*(1.0 +77.0),.001);
//		
//	}
	
	// Test potassium time independent current, IK1
	void dontTestPotassiumTimeIndependentCurrentLR91(void)
	{
		PotassiumTimeIndependentCurrentLR91 *myPotassiumTimeIndependentCurrent;		
		myPotassiumTimeIndependentCurrent = new PotassiumTimeIndependentCurrentLR91();
		double v = 2.0 ;
		myPotassiumTimeIndependentCurrent->UpdateAlphaAndBeta(v) ;
		myPotassiumTimeIndependentCurrent->UpdateMagnitudeOfCurrent(v) ;
		double K1test = myPotassiumTimeIndependentCurrent->GetK1(v) ;
		TS_ASSERT_DELTA(myPotassiumTimeIndependentCurrent->GetMagnitudeOfCurrent(), .6047*K1test*(v- (-77)) , .001) ;
		
		
		//TS_ASSERT(true);
	}
	
	// Test plateau potassium current, IKp
	void dontTestPlateauPotassiumCurrentLR91(void)
	{
		PlateauPotassiumCurrentLR91 *myPlateauPotassiumCurrent;
		myPlateauPotassiumCurrent = new PlateauPotassiumCurrentLR91();
		TS_ASSERT(true);
	}
	
	// Test background current, IB
	void dontTestBackgroundCurrentLR91(void)
	{
		// std::cout << "Running TestBackgroundCurrentLR91..." << std::endl;
		BackgroundCurrentLR91 *myBackgroundCurrent;
		myBackgroundCurrent = new BackgroundCurrentLR91();
		TS_ASSERT(true);
	}
};

#endif //_TESTLR91CURRENTS_HPP_
