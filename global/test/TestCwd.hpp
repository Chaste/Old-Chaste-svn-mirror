#ifndef _TESTCWD_HPP_
#define _TESTCWD_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdlib>

#include "PetscSetupAndFinalize.hpp"

/**
 * Test for a strange 'feature' of Debian sarge systems, where the
 * current working directory changes on PETSc initialisation.
 *
 * Define CWD_HACK to work around this (see finarfin and maths systems
 * for examples).
 */
class TestCwd : public CxxTest::TestSuite
{
public:
    void TestShowCwd()
    {
        TS_ASSERT_EQUALS(system("pwd"), 0);
        TS_ASSERT_EQUALS(system("ls -l io/test/data"), 0);
    }
};

#endif
