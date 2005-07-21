#ifndef _PETSCSETUPANDFINALIZE_HPP_
#define _PETSCSETUPANDFINALIZE_HPP_

/**
 * This file is designed to be included by any test suites that use PETSc.
 * It controls the PETSc initialisation and finalisation.
 * 
 * Currently it will dump info on any non-freed vectors or matrices
 * on finalisation.
 */

#include <cxxtest/GlobalFixture.h>
#include "petsc.h"  

class PetSCSetup : public CxxTest::GlobalFixture
{
public:
	/// Standard setup method for PETSc.
	bool setUpWorld()
	{
#ifdef PETSC_MEMORY_TRACING
		int FakeArgc = 4;
		char *FakeArgs[] = {"testrunner", "-trmalloc", "-trdebug", "-trdump"};
#else
		int FakeArgc = 1;
		char *FakeArgs[] = {"testrunner"};
#endif
		char **FakeArgv=(char **)FakeArgs;

		PetscErrorCode ierr = PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, PETSC_NULL);
		
		//std::cout << std::endl << "Petsc init: " << ierr << std::endl;
		
		return true;
	}
	// Clean up PETSc after running all tests.
	bool tearDownWorld() 
    {
    	//std::cout << "Finalizing..." << std::endl;
    	PetscFinalize();
    	return true;
    }
};
static PetSCSetup thisSetup;


#endif //_PETSCSETUPANDFINALIZE_HPP_
