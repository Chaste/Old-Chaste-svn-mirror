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
		 //std::cout << "\n Running TestPotassiumTimeDependentCurrentLR91...\n \n" << std::endl;
		PotassiumTimeDependentCurrentLR91 *myPotassiumTimeDependentCurrent;
		myPotassiumTimeDependentCurrent = new PotassiumTimeDependentCurrentLR91();
		
		myPotassiumTimeDependentCurrent->SetGatingVariables(1.0);
		myPotassiumTimeDependentCurrent->UpdateMagnitudeOfCurrent(1.0, 1.0);
		double xtest=myPotassiumTimeDependentCurrent->GetX() ;
		double xitest =myPotassiumTimeDependentCurrent->GetXi(1.0) ;
		TS_ASSERT(myPotassiumTimeDependentCurrent->GetX() == 1.0);
		//TS_ASSERT(myPotassiumTimeDependentCurrent->GetMagnitudeOfCurrent() == 1.0);
		std::cout<<"MAgnitude of current is :" << myPotassiumTimeDependentCurrent->GetMagnitudeOfCurrent() << "\n \n ...."  ;
		std::cout<<gK << " \n " << xtest << " \n " <<xitest << " \n " ;
		//TS_ASSERT_DELTA(myPotassiumTimeDependentCurrent->GetMagnitudeOfCurrent() , gK*xtest*xitest*(1.0 +77.0),.02);
		
	}
	
	// Test potassium time independent current, IK1
	void TestPotassiumTimeIndependentCurrentLR91(void)
	{
		PotassiumTimeIndependentCurrentLR91 *myPotassiumTimeIndependentCurrent;
		
		myPotassiumTimeIndependentCurrent = new PotassiumTimeIndependentCurrentLR91();
		
		myPotassiumTimeIndependentCurrent->UpdateAlphaAndBeta(2.0) ;
		myPotassiumTimeIndependentCurrent->UpdateMagnitudeOfCurrent(2.0) ;
		double I = myPotassiumTimeIndependentCurrent->GetMagnitudeOfCurrent() ;
		double v = 2.0 ;
		double K1test = myPotassiumTimeIndependentCurrent->GetK1(2.0) ;
		
		std::cout<< " The value of K1test is............. \n"<<  K1test  << "\n \n \n " ;
		TS_ASSERT_DELTA(myPotassiumTimeIndependentCurrent->GetMagnitudeOfCurrent(), .6047*K1test*(v- (-77)) , .001) ;
		
		
		//TS_ASSERT(true);
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
