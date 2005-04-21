#ifndef _TESTABSTRACTODESYSTEM_HPP_
#define _TESTABSTRACTODESYSTEM_HPP_

// TestAbstractOdeSystem.hpp

#include <cmath>
#include "petscvec.h"
#include "petsc.h"
#include <iostream>
#include <vector>
#include "AbstractOdeSystem.hpp"
#include "TestOde1.hpp"
#include "TestOde2.hpp"
#include "TestOde3.hpp"

// Tolerance for tests
double tol=0.01;

class TestAbstractOdeSystem : public CxxTest::TestSuite
{
	public:
	
	void setUp()
	{
		PetscInitialize(0,NULL,0,0);
	}
	
	void tearDown()
	{
	}
		
	void testAddition(void)
	{
		std::cout << "\n Start of OdeSystem tests\n";
		
		TS_ASSERT(1+1 >1)
	}
	
	void testPetsc(void)
	{
		int SystemSize=4;
		PetscScalar Value = 2.5;
		PetscScalar * yInitElements;
		std::cout <<  SystemSize << "\n";
		Vec  yInit;
		
		// Creation of vector with Petsc
		VecCreate(PETSC_COMM_WORLD,&yInit);
		VecSetSizes(yInit,PETSC_DECIDE,SystemSize);
		VecSetType(yInit,VECSEQ);
		// Set all values in the vector to Value
		VecSet(&Value, yInit);
		// Convert PetscVector to pointer
		VecGetArray(yInit, &yInitElements);
		for(int k=0; k < SystemSize; k++)
		{
			std::cout << yInitElements[k] << "\t";
		}
		std::cout << "\n";
		yInitElements[2]=5.0;
		TS_ASSERT(yInitElements[2] > 4);
		VecRestoreArray(yInit,&yInitElements);
		// Reading(Viewing) in a PetScvector
		VecView(yInit,PETSC_VIEWER_STDOUT_WORLD);
		
		//std::cout << "\n End of first set of tests\n";
	}
	
	void TestMyAbstractOdeSystem(void)
	{
		//std::cout << "Start of second test\n";
		int nEqus=1;
		double * yInit = new double[nEqus];
		double value = 1.0;
		
		for(int k=0; k<nEqus; k++)
		{
			yInit[k]= value;
		}
		
		double StartTime= 0.0;
		
		//commented out as EvaluateYPrime is defined in subclasses and ruins test
/*		AbstractOdeSystem Mine(nEqus,StartTime, yInit);
		std::cout << "\n";
		for(int k=0; k<nEqus; k++)
		{
			std::cout << yInit[k] << "\n";
		}
		// Finish the test
*/		//std::cout << "\n End of second set of tests\n";
	}
	
	void TestOdeOne(void)
	{
		std::cout << "Start of third test\n";
		std::vector<double> yInit(1);
		yInit[0]=0.0;
		
//		TestOde1 Ode1();
//		double rYPrime[1];
//		Ode1.EvaluateYDerivatives(1.0, yInit, rYPrime);
//		std::cout << rYPrime[0];
		
		TestOde1 * pOde1 = new TestOde1();
		std::vector<double> rYPrime2;
		rYPrime2=pOde1->EvaluateYDerivatives(1.0, yInit);
		std::cout << rYPrime2[0];
		TS_ASSERT_DELTA(rYPrime2[0],1.0,tol);
		std::cout << "\n End of third set of tests\n";
	}
	
	
	void TestOdeTwo(void)
	{
		std::cout << "Start of fourth test\n";
		std::vector<double> yInit(1);
		
		yInit[0]=4.0;
		
		TestOde2 * pOde2 = new TestOde2();
		std::vector<double> rYPrime2;
		rYPrime2=pOde2->EvaluateYDerivatives(2.0, yInit);
		std::cout << rYPrime2[0];
		TS_ASSERT_DELTA(rYPrime2[0],8.0,tol);
		std::cout << "\n End of fourth set of tests\n";
	}
	
	void TestOdeThree(void)
	{
		std::cout << "Start of fifth test\n";
		std::vector<double> yInit(2);
		
		yInit[0]=4.0;
		yInit[1]=8.0;
		
		TestOde3 * pOde3 = new TestOde3();
		std::vector<double> rYPrime3;
		rYPrime3=pOde3->EvaluateYDerivatives(2.0, yInit);
		std::cout << rYPrime3[0] << std::endl;
		std::cout << rYPrime3[1] << std::endl;
		TS_ASSERT_DELTA(rYPrime3[0],8.0,tol);
		TS_ASSERT_DELTA(rYPrime3[1],16.0,tol);
		std::cout << "\n End of fifth set of tests\n";
		
		
	std::cout << "\n End of OdeSystem tests\n";
	}
	
};



#endif //_TESTABSTRACTODESYSTEM_HPP_
