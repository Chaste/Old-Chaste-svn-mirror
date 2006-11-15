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

#ifdef CWD_HACK
#  include <unistd.h>
#  include <iostream>
#endif

#include "PetscException.hpp"

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
        PETSCEXCEPT(PetscInitialize(&cxxtest_argc, &cxxtest_argv,
                                    PETSC_NULL, PETSC_NULL) );
                                    
#ifdef CWD_HACK
        char buf[10000];
        std::cout << std::endl << "CWD: " << getcwd(buf, 10000) << std::endl;
        chdir("../../..");
        std::cout << "CWD: " << getcwd(buf, 10000) << std::endl;
#endif
                                    
        return true;
    }
    // Clean up PETSc after running all tests.
    bool tearDownWorld()
    {
        //std::cout << "Finalizing..." << std::endl;
        PETSCEXCEPT(PetscFinalize());
        return true;
    }
};
static PetSCSetup thisSetup;


#endif //_PETSCSETUPANDFINALIZE_HPP_
