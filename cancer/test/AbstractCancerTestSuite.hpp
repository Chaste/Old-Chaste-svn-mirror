#ifndef ABSTRACTCANCERTESTSUITE_HPP_
#define ABSTRACTCANCERTESTSUITE_HPP_

/**
 * This class provides setUp and tearDown methods that are common to
 * many cancer test suites.  Such suites may inherit from this class
 * to avoid having to redefine them.
 */
class AbstractCancerTestSuite : public CxxTest::TestSuite
{
protected:    
    void setUp()
    {
        // Initialise singleton classes
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);
        CancerParameters::Instance()->Reset();
    }
    
    void tearDown()
    {
        // Clear up singleton classes
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};


#endif /*ABSTRACTCANCERTESTSUITE_HPP_*/
