#ifndef RANDOMNUMBERGENERATORS_HPP_
#define RANDOMNUMBERGENERATORS_HPP_
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>

class RandomNumberGenerator
{
public:
    double StandardNormalRandomDeviate(void);
    double NormalRandomDeviate(double mean, double sd);
    double ranf(void);
    unsigned randMod(unsigned base);
        
    static RandomNumberGenerator* Instance();
    static void Destroy();
    void Reseed(int seed)
    {
        mSeed = seed;
        srandom(seed);
        mTimesCalled = 0;
    }
protected:
    /**
      * @param seed Is the new seed which defaults to zero.
      */
    RandomNumberGenerator()
    {
        mSeed = 0;
        mTimesCalled = 0;
        srandom(0);
    }
    
private:
    int mSeed;
    unsigned mTimesCalled;

    static RandomNumberGenerator* mpInstance;
    
    friend class boost::serialization::access;
    /**
     * Serialization of a RandomNumberGenerator object must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        // note, version is always the latest when saving
        archive & mSeed;
        archive & mTimesCalled;
    }
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & mSeed;
        archive & mTimesCalled;
        // reset the random number generator to use the correct seed
        srandom(mSeed);
        // call it the correct number of times to put it in the 
        // same state as it used to be.
        // NOTE: This is only guaranteed to work if Normal random
        // deviates are not used, since the methods to generate
        // numbers from a normal distribution use static variables.
        for (unsigned i=0; i<mTimesCalled; i++)
        {
            random();
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    
};
#endif /*RANDOMNUMBERGENERATORS_HPP_*/
