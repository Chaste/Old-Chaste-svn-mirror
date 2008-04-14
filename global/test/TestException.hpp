/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
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
        
        TS_ASSERT_THROWS_ANYTHING(EXCEPTION("Hello. I'm an exception"));
        TS_ASSERT_THROWS_ANYTHING(NEVER_REACHED);
    }
};

#endif //_TESTEXCEPTION_HPP_
