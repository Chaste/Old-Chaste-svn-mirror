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
    }
};

#endif //_TESTEXCEPTION_HPP_
