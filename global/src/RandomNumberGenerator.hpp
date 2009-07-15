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

#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>

/**
 * A special singleton class allowing one to generate different types of
 * random number in a globally consistent way.
 */
class RandomNumberGenerator
{
public:


    /**
     * Generate a random number from the normal distribution with mean 0
     * and standard distribution 1.
     */
    double StandardNormalRandomDeviate();

    /**
     * Generate a random number from a normal distribution with given
     * mean and standard deviation.
     *
     * @param mean the mean of the normal distribution from which the random number is drawn
     * @param sd the standard deviation of the normal distribution from which the random number is drawn
     */
    double NormalRandomDeviate(double mean, double sd);

    /**
     * Generate a uniform random number in (0,1).
     */
    double ranf();

    /**
     * Generate a random number modulo base (ie an integer
     * within the range 0,..,base-1),
     *
     * @param base the order of the field of positive integers from which the random number is drawn
     */
    unsigned randMod(unsigned base);

    /**
     * Shuffle the integers 0,1,..,num-1, using the Knuth-algorithm
     * (also called the Fisher-Yates algorithm), a linear time unbiased method.
     * The shuffled values are returned in rValues, which doesn't need to
     * be correctly-sized when passed in.
     *
     * @param num  the number of integers to shuffle
     * @param rValues  the shuffled values
     */
    void Shuffle(unsigned num, std::vector<unsigned>& rValues);

    /**
     * Return a pointer to the random number generator object.
     * The object is created the first time this method is called.
     */
    static RandomNumberGenerator* Instance();

    /**
     * Destroy the current instance of the random number generator.
     * The next call to Instance will create a new instance and re-seed.
     * This method *must* be called before program exit, to avoid a memory
     * leak.
     */
    static void Destroy();

    /**
     * Reseed the random number generator.
     *
     * @param seed the new seed
     */
    void Reseed(int seed);

protected:

    /**
     * Protected constructor.
     * Use Instance() to access the random number generator.
     */
    RandomNumberGenerator();

private:

    /** The random number generator seed. */
    int mSeed;

    /** The number of times the random number generator has been called. */
    unsigned mTimesCalled;

    /** Pointer to the single instance. */
    static RandomNumberGenerator* mpInstance;

    friend class boost::serialization::access;
    /**
     * Save the RandomNumberGenerator and its member variables.
     *
     * Serialization of a RandomNumberGenerator object must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        // note, version is always the latest when saving
        archive & mSeed;
        archive & mTimesCalled;
    }
    /**
     * Load the RandomNumberGenerator and its member variables.
     *
     * Serialization of a RandomNumberGenerator object must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
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
