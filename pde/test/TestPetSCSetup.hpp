#ifndef _TESTPETSCSETUP_HPP_
#define _TESTPETSCSETUP_HPP_

#include <cxxtest/TestSuite.h>
#include <cxxtest/GlobalFixture.h>
#include <petscvec.h>
#include <petscmat.h>


/*
 * 
 * This test was useful when we were using many test header files per test executable.
 * The idea was to make sure that the first invocation of  PetscInitialize had
 * the correct arguments.  It should now be unnecessary, but is included for completeness.
 * See below for alternative setUpWorld etc.
 */
 
 
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

#ifdef __THIS_IS_AN_ALTERNATIVE_CONFIGURATION__

	bool setUpWorld(){ 
  		/*
  		 */
  		int FakeArgc=1;
		char *FakeArgs[3];
		char *FakeArgv0="testrunner";
		FakeArgs[0]=FakeArgv0;
		char *FakeArgv1="-trmalloc_log";
		FakeArgs[1]=FakeArgv1;
		char *FakeArgv2="-trdebug";
		FakeArgs[2]=FakeArgv2;
		
		
        
        char **FakeArgv=(char **)FakeArgs;
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
    	/*
    	*/
    	
    	
    	//PetscInitializeNoArguments();
 			/*
 			 * 
 			PetscLogAllBegin();
         	PetscTrLog();
        *
        */
    	return true; 
 	}
    bool tearDownWorld() 
    { 
/*    FILE *p;
    p=fopen("temp","w");
    	PetscTrLogDump(p);
    	*
    	*/
    	PetscFinalize(); 
    	return true;
    }
#endif
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
