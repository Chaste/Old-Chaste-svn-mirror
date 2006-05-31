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
#include <petsc.h>
#include <stdlib.h>


#define PETSCEXCEPT(n) if (n) PetscException(n, __LINE__, __FUNCT__,__FILE__);

void PetscException(int petscError,
                    int line,
                    const char* funct, 
                    const char* file)
{
    const char*  text;
 
    PetscErrorMessage(petscError,  &text, NULL);
    
    std::string err_string;
    err_string=text;
    err_string+= " in function ";
    err_string+= funct;
    err_string+= " on line " ;
    err_string+= line;
    err_string+= " of file ";
    err_string+= file;
    
    Exception e(err_string);
    throw(e);
}

class PetSCSetup : public CxxTest::GlobalFixture
{
public:
	/// Standard setup method for PETSc.
	bool setUpWorld()
	{
        /**
         * The cxxtest_argc and cxxtest_argv variables are global, and filled in
         * from the arguments passed to the cxxtest test suite runner.
         */
		//PetscErrorCode ierr = 
        PetscInitialize(&cxxtest_argc, &cxxtest_argv, PETSC_NULL, PETSC_NULL);
		
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
