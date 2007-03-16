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
        
        RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        
        double ran2=p_gen->ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
        
        RandomNumberGenerator::Destroy();
    }
    
    void TestNewMethodSeed()
    {
        RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        double ran2=p_gen->ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
        
        RandomNumberGenerator::Destroy();
        
    }
    
    void TestDifferentRandomSeed()
    {
        srandom(36);
        ran1=(double)random()/RAND_MAX;
        
        RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        p_gen->Reseed(36);
        
         
        double ran2=p_gen->ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);        

        RandomNumberGenerator::Destroy();

    }
};

#endif /*TESTRANDOMNUMBERS_HPP_*/
