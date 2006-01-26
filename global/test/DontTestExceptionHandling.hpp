#ifndef _TESTEXCEPTIONHANDLING_HPP_
#define _TESTEXCEPTIONHANDLING_HPP_

#include <cxxtest/TestSuite.h>
#include "Exception.hpp"

/**
 * This tests that cxxtest handles exceptions thrown in tests gracefully.
 */
 
 
class TestPetSCSetup : public CxxTest::TestSuite
{
public:
    void TestThrowingAnExceptionInATest()
    {
        throw Exception("Will cxxtest be nice I wonder?");
    }
};

#endif //_TESTEXCEPTIONHANDLING_HPP_
