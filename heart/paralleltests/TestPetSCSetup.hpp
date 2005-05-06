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
