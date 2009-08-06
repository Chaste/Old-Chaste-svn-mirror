/*

Copyright (C) University of Oxford, 2005-2009

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

#ifndef _TESTEXCEPTION_HPP_
#define _TESTEXCEPTION_HPP_

#include <cxxtest/TestSuite.h>
#include "Exception.hpp"

class TestException : public CxxTest::TestSuite
{
public:
    void TestGetMessage()
    {
        std::string msg("This is an exception");

        try
        {
            EXCEPTION(msg);
        }
        catch (Exception e)
        {
            std::string e_msg = e.GetMessage();
            std::string::size_type e_len = e_msg.length();
            std::string::size_type len = msg.length();
            TS_ASSERT_EQUALS(e_msg.substr(e_len - len), msg);
        }

        TS_ASSERT_THROWS_EQUALS(EXCEPTION("Hello. I'm an exception"), const Exception &err,
                        err.GetShortMessage(), "Hello. I'm an exception" );
        //NB The following test will fail if the number of lines above is changed... (that's why method GetShortMessage() was introduced).
        TS_ASSERT_THROWS_EQUALS(EXCEPTION("Hello. I'm an exception"), const Exception &err,
                                err.GetMessage(), "\nChaste error: ./global/test/TestException.hpp:58: Hello. I\'m an exception" );

        TS_ASSERT_THROWS_EQUALS(NEVER_REACHED,  const Exception &err,
                err.GetShortMessage(), "Should have been impossible to reach this line of code");
    }
};

#endif //_TESTEXCEPTION_HPP_
