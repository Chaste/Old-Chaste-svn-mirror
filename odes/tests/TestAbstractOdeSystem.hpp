#ifndef _TESTABSTRACTODESYSTEM_HPP_
#define _TESTABSTRACTODESYSTEM_HPP_

// TestAbstractOdeSystem.hpp

#include <cmath>
//#include "petscvec.h"
//#include "petsc.h"
#include <iostream>
#include "AbstractOdeSystem.hpp"
#include "TestOde1.hpp"

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
		
		std::cout << "\n End of first set of tests\n";
	}
	
	void TestMyAbstractOdeSystem(void)
	{
		std::cout << "Start of second test\n";
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
*/		std::cout << "\n End of second set of tests\n";
	}
	
	void TestOdeOne(void)
	{
		double* yInit = new double[1];
		yInit[0]=0.0;
		
		TestOde1 Ode1(0.0, yInit);
		double rYPrime[1];
		Ode1.EvaluateYPrime(1.0, yInit, rYPrime);
		std::cout << rYPrime[0];
		
		TestOde1 * pOde1 = new TestOde1(0.0, yInit);
		double rYPrime2[1];
		pOde1->EvaluateYPrime(1.0, yInit, rYPrime2);
		std::cout << rYPrime2[0];
	}
};



#endif //_TESTABSTRACTODESYSTEM_HPP_
