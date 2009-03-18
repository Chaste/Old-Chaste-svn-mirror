#ifndef TESTHELLO_HPP_
#define TESTHELLO_HPP_

#include <cxxtest/TestSuite.h>

#include "Hello.hpp"

class TestHello : public CxxTest::TestSuite
{
public:
    void testHello()
    {
        Hello world("Hello world!");
        TS_ASSERT_EQUALS(world.GetMessage(), "Hello world!");
        TS_ASSERT_THROWS_ANYTHING(world.Complain("I don't like you"));
    }
};

#endif /*TESTHELLO_HPP_*/
