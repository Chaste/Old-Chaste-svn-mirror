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


#ifndef _CHASTEPOINT_HPP_
#define _CHASTEPOINT_HPP_

#include "UblasIncludes.hpp"
#include <cassert>
#include <vector>
#include "Exception.hpp"
using namespace boost::numeric::ublas;

/**
 * A ChastePoint class, templated over spatial dimension.
 */
template<unsigned DIM>
class ChastePoint
{
private:

    /** The location of the Point. */
    c_vector<double, DIM> mLocation;

public:

    /**
     * Create a Point object.
     * There are 3 optional arguments, which can be used to specify the values
     * of the first 3 dimensions, if present.
     *
     * Point now uses a ublas vector to store its location. The
     * rGetLocation method returns a reference to this vector.
     * Use of this method together with ublas operations
     * is the perfered way to use this class.
     * 
     * @param v1  the point's x-coordinate (defaults to 0)
     * @param v2  the point's y-coordinate (defaults to 0)
     * @param v3  the point's z-coordinate (defaults to 0)
     */
    ChastePoint(double v1=0, double v2=0, double v3=0)
    {
        if (DIM>0)
        {
            mLocation[0] = v1;
        }
        if (DIM>1)
        {
            mLocation[1] = v2;
        }
        if (DIM>2)
        {
            mLocation[2] = v3;
        }
    }

    /**
     * Create a Point object.
     * This constructor takes a vector giving the coordinates of the point.
     * The length of the vector must be at least the dimension of the point.
     * 
     * @param coords  a std::vector storing the point's coordinates
     */
    ChastePoint(std::vector<double> coords)
    {
        for (unsigned i=0; i<DIM; i++)
        {
            mLocation(i) = coords.at(i);
        }
    }

    /**
     * Alternative constructor which takes in a c_vector.
     * 
     * @param location  a c_vector storing the point's coordinates
     */
    ChastePoint(c_vector<double, DIM> location)
    {
        mLocation = location;
    }

    /**
     * Get the location of the Point.
     */
    c_vector<double, DIM>& rGetLocation(void)
    {
        return mLocation;
    }

    /**
     * Get the location of the Point.
     */
    const c_vector<double, DIM>& rGetLocation(void) const
    {
        return mLocation;
    }

    /**
     * Access the vector mLocation.
     * 
     * @param i the index of the vector to return
     */
    double operator[] (unsigned i) const
    {
        assert(i<DIM);
        return mLocation(i);
    }

    /**
     * Set one of the coordinates of the Point.
     * 
     * @param i the index of the coordinate
     * @param value the value of the coordinate
     */
    void SetCoordinate(unsigned i, double value)
    {
        assert(i<DIM);
        mLocation(i) = value;
    }
};

/**
 * A  zero-dimensional ChastePoint class.
 */
template<>
class ChastePoint<0>
{
public:

    /**
     * Create a zero-dimensional Point object.
     * There are 3 optional arguments, which should not be used.
     * 
     * @param v1  the point's x-coordinate (defaults to 0)
     * @param v2  the point's y-coordinate (defaults to 0)
     * @param v3  the point's z-coordinate (defaults to 0)
     */
    ChastePoint(double v1=0, double v2=0, double v3=0)
    {
    }

    /**
     * Access the vector mLocation.
     * 
     * @param i the index of the vector to return
     */
    double operator[] (unsigned i) const
    {
        EXCEPTION("Zero-dimensional point has no data");
    }
};


#endif //_CHASTEPOINT_HPP_
