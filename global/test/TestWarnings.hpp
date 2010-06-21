/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef _TESTWARNINGS_HPP_
#define _TESTWARNINGS_HPP_

#include <cxxtest/TestSuite.h>
#include "Warnings.hpp"

class TestWarnings: public CxxTest::TestSuite
{
public:
    void TestGetMessage()
    {
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        WARNING("Ozzy is near.");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        WARNING("Stay alert.");   
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Ozzy is near.");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        Warnings::QuietDestroy();
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
    }
    
    void TestWarningsNoisy()
    {
        TS_ASSERT_THROWS_THIS(Warnings::Instance()->GetNextWarningMessage(),"There are no warnings");
        WARNING("This one will get printed.");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
    }
        

 };

#endif //_TESTWARNINGS_HPP_
