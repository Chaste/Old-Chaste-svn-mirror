#ifndef _TESTPETSCSETUP_HPP_
#define _TESTPETSCSETUP_HPP_

#include <cxxtest/TestSuite.h>
#include <cxxtest/GlobalFixture.h>
#include <petscvec.h>
#include <petscmat.h>


class PetSCSetup : public CxxTest::GlobalFixture
{
public:
  
 	bool setUpWorld(){ 
  		/*
  		 * 
  		int FakeArgc=0;
		char *FakeArgv0="testrunner";
		char **FakeArgv=&FakeArgv0;
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
    	*
    	*/
    	
    	
    	PetscInitializeNoArguments();
        
    	return true; 
 	}
    bool tearDownWorld() 
    { 
    	PetscFinalize(); 
    	return true;
    }
};

static PetSCSetup thisSetup;

class TestPetSCSetup : public CxxTest::TestSuite
{
public:
    void testPetscIsThere()
    {
        PetscTruth is_there;
        PetscInitialized(&is_there);
        TS_ASSERT( is_there == PETSC_TRUE );
    }
};



#endif // _TESTPETSCSETUP_HPP_
