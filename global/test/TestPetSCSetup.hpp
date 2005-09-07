#ifndef _TESTPETSCSETUP_HPP_
#define _TESTPETSCSETUP_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"

/**
 * This tests that the initialisation of PETSc does something.
 */
 
 
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
