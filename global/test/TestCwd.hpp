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
