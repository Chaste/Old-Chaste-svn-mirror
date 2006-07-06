#ifndef _TESTEXCEPTIONHANDLING_HPP_
#define _TESTEXCEPTIONHANDLING_HPP_

#include <cxxtest/TestSuite.h>
#include "Exception.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * This tests that cxxtest handles exceptions thrown in tests gracefully.
 *
 * The tests are supposed to fail, so aren't run normally.
 */
 
 
class TestExceptionHandling : public CxxTest::TestSuite
{
public:
    void TestThrowingWithPetsc()
    {
        Vec test_vec;
        VecCreate(PETSC_COMM_WORLD, &test_vec);
        VecSetSizes(test_vec, PETSC_DECIDE, 20);
        VecSetFromOptions(test_vec);
        
        VecAssemblyBegin(test_vec);
        VecAssemblyEnd(test_vec);
        
        throw Exception("Will cxxtest be nice if we do PETSc things?");
    }
    
    void ThrowExceptionMethod()
    {
        throw Exception("Exception thrown from method.");
    }

    void TestThrowingAnExceptionInATest()
    {
        throw Exception("Will cxxtest be nice I wonder?");
    }

    void TestCatchingExceptionWithCxxtest()
    {
        TS_ASSERT_THROWS_ANYTHING(throw Exception("Will cxxtest be nice I wonder?"));
        TS_ASSERT_THROWS_NOTHING(throw Exception("Will cxxtest be nice I wonder?"));
    }
};

#endif //_TESTEXCEPTIONHANDLING_HPP_
