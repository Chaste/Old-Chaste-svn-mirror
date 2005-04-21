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
	
	// Test potassium time dependent current, IK
	void TestPotassiumTimeDependentCurrentLR91(void)
	{
		// std::cout << "Running TestPotassiumTimeDependentCurrentLR91..." << std::endl;
		PotassiumTimeDependentCurrentLR91 *myPotassiumTimeDependentCurrent;
		
		myPotassiumTimeDependentCurrent = new PotassiumTimeDependentCurrentLR91();
		myPotassiumTimeDependentCurrent->SetGatingVariables(1.0);
		myPotassiumTimeDependentCurrent->SetMagnitudeOfCurrent(50.0);
		
		TS_ASSERT(myPotassiumTimeDependentCurrent->GetX() == 1.0);
		TS_ASSERT(myPotassiumTimeDependentCurrent->GetMagnitudeOfCurrent() == 50.0);
	}
	
	// Test potassium time independent current, IK1
	void TestPotassiumTimeIndependentCurrentLR91(void)
	{
		PotassiumTimeIndependentCurrentLR91 *myPotassiumTimeIndependentCurrent;
		
		myPotassiumTimeIndependentCurrent = new PotassiumTimeIndependentCurrentLR91();
		
		TS_ASSERT(true);
	}
	
	// Test plateau potassium current, IKp
	void TestPlateauPotassiumCurrentLR91(void)
	{
		PlateauPotassiumCurrentLR91 *myPlateauPotassiumCurrent;
		
		myPlateauPotassiumCurrent = new PlateauPotassiumCurrentLR91();
		
		TS_ASSERT(true);
	}
	
	// Test background current, IB
	void TestBackgroundCurrentLR91(void)
	{
		// std::cout << "Running TestBackgroundCurrentLR91..." << std::endl;
		BackgroundCurrentLR91 *myBackgroundCurrent;
		
		myBackgroundCurrent = new BackgroundCurrentLR91();
		
		TS_ASSERT(true);
	}
};

#endif //_TESTLR91CURRENTS_HPP_
