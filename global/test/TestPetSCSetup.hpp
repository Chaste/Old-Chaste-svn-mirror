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
    void TestPetscIsThere()
    {
        PetscTruth is_there;
        PetscInitialized(&is_there);
        TS_ASSERT( is_there == PETSC_TRUE );
    }
    
    void TestPetscExceptions()
    {
        int err=0;
        TS_ASSERT_THROWS_NOTHING(PETSCEXCEPT(err));
        
        Vec v;
        err = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, -1, &v);
        VecDestroy(v);
        //#define PETSC_ERR_ARG_WRONGSTATE   73   /* object in argument is in wrong */
        TS_ASSERT_THROWS_ANYTHING(PETSCEXCEPT(err));
        
        err=PETSC_ERR_FILE_OPEN;
        //#define PETSC_ERR_FILE_OPEN        65   /* unable to open file */
        TS_ASSERT_THROWS_ANYTHING(PETSCEXCEPT(err));
        
        //See if we can do it without a temporary
        TS_ASSERT_THROWS_ANYTHING(
            PETSCEXCEPT(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, -1, &v)));
        VecDestroy(v);
            
        //This test give back an "unknown error" message
        TS_ASSERT_THROWS_ANYTHING( PETSCEXCEPT(-3));
    }
    
    
    void TestKspExceptionsForCoverage()
    {
        TS_ASSERT_THROWS_NOTHING(  KSPEXCEPT(2) );
        //These next few lines are designed to force the coverage test to pass.
        //Some are hard to throw in normal circumstances --
        //"Unknown KSP error code" ought never to be thrown.
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_ITS) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_DTOL) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_BREAKDOWN) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_BREAKDOWN_BICG) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_NONSYMMETRIC) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_INDEFINITE_PC) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(-735827) );
    }
};



#endif // _TESTPETSCSETUP_HPP_
