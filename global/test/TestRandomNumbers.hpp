#ifndef TESTRANDOMNUMBERS_HPP_
#define TESTRANDOMNUMBERS_HPP_
#include <cxxtest/TestSuite.h>

#include "RandomNumberGenerators.hpp"

class TestRandomNumbers : public CxxTest::TestSuite
{
public:
    double ran1;
    void TestRandomNumers()
    {
        srandom(0);
        ran1=(double)random()/RAND_MAX;
        
        RandomNumberGenerators gen;
        
        double ran2=gen.ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
        
        
    }
    void TestNewMethodSeed()
    {
        RandomNumberGenerators gen;
        double ran2=gen.ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
    }
    void TestDifferentRandomSeed()
    {
        srandom(36);
        ran1=(double)random()/RAND_MAX;
        
        RandomNumberGenerators gen(36);
        
        double ran2=gen.ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
        
        
    }
};
#endif /*TESTRANDOMNUMBERS_HPP_*/
