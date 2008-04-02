/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

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
        
        EXCEPTION("Will cxxtest be nice if we do PETSc things?");
    }
    
    void ThrowExceptionMethod()
    {
        EXCEPTION("Exception thrown from method.");
    }
    
    void TestThrowingAnExceptionInATest()
    {
        EXCEPTION("Will cxxtest be nice I wonder?");
    }
    
    void TestCatchingExceptionWithCxxtest()
    {
        TS_ASSERT_THROWS_ANYTHING(EXCEPTION("Will cxxtest be nice I wonder?"));
        TS_ASSERT_THROWS_NOTHING(EXCEPTION("Will cxxtest be nice I wonder?"));
    }
};

#endif //_TESTEXCEPTIONHANDLING_HPP_
