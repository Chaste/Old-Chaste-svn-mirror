#ifndef _TESTPETSCSETUP_HPP_
#define _TESTPETSCSETUP_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

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
    
    void testPetscExceptions()
    {
        int err=0;
        TS_ASSERT_THROWS_NOTHING(PETSCEXCEPT(err));
       
        Vec v;     
        err = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, -1, &v);
        //#define PETSC_ERR_ARG_WRONGSTATE   73   /* object in argument is in wrong */
        TS_ASSERT_THROWS_ANYTHING(PETSCEXCEPT(err));
        
        err=PETSC_ERR_FILE_OPEN;
        //#define PETSC_ERR_FILE_OPEN        65   /* unable to open file */
        TS_ASSERT_THROWS_ANYTHING(PETSCEXCEPT(err));
        
        //See if we can do it without a temporary
        TS_ASSERT_THROWS_ANYTHING(
            PETSCEXCEPT(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, -1, &v)));
        
        //This test give back an "unknown error" message
        TS_ASSERT_THROWS_ANYTHING( PETSCEXCEPT(-3));    
    }
};



#endif // _TESTPETSCSETUP_HPP_
