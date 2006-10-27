#ifndef TESTRANDOMNUMBERS_HPP_
#define TESTRANDOMNUMBERS_HPP_
#include <cxxtest/TestSuite.h>

#include "RandomNumberGenerator.hpp"

class TestRandomNumbers : public CxxTest::TestSuite
{
public:
    double ran1;

    void TestRandomNumers()
    {
        srandom(0);
        ran1=(double)random()/RAND_MAX;
        
        RandomNumberGenerator gen;
        
        double ran2=gen.ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
    }
    
    void TestNewMethodSeed()
    {
        RandomNumberGenerator gen;
        double ran2=gen.ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
    }
    
    void TestDifferentRandomSeed()
    {
        srandom(36);
        ran1=(double)random()/RAND_MAX;
        
        RandomNumberGenerator gen(36);
        
        double ran2=gen.ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);        
    }
};

#endif /*TESTRANDOMNUMBERS_HPP_*/
