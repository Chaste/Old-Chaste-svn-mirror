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


#ifndef RANDOMNUMBERGENERATORS_HPP_
#define RANDOMNUMBERGENERATORS_HPP_
#include <cmath>
#include <ctime>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>

/**
 * A special singleton class allowing one to generate different types of 
 * random number in a globally consistent way.
 */
class RandomNumberGenerator
{
public:
    double StandardNormalRandomDeviate();
    double NormalRandomDeviate(double mean, double sd);
    double ranf();
    unsigned randMod(unsigned base);
    void Shuffle(unsigned num, std::vector<unsigned>& rValues);


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
