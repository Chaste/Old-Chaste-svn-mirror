/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _CHASTEPOINT_HPP_
#define _CHASTEPOINT_HPP_


#include <boost/numeric/ublas/vector.hpp>
#include <cassert>
#include <vector>
#include "Exception.hpp"
using namespace boost::numeric::ublas;


template<unsigned DIM>
class ChastePoint
{
private:
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
     */
    ChastePoint(std::vector<double> coords)
    {
        for (unsigned i=0; i<DIM; i++)
        {
            mLocation(i) = coords.at(i);
        }
    }
    
    ChastePoint(c_vector<double, DIM> location)
    {
        mLocation = location;
    }
    
    c_vector<double, DIM>& rGetLocation(void)
    {
        return mLocation;
    }
    
    const c_vector<double, DIM>& rGetLocation(void) const
    {
        return mLocation;
    }
    
    double operator[] (unsigned i) const
    {
        assert(i<DIM);
        return mLocation(i);
    }
    
    void SetCoordinate(unsigned i, double value)
    {
        assert(i<DIM);
        mLocation(i) = value;
    }
};

template<>
class ChastePoint<0>
{
public:

    /**
     * Create a zero-dimensional Point object.
     * There are 3 optional arguments, which should not be used.
     */
    ChastePoint(double v1=0, double v2=0, double v3=0)
    {
    }
    
    double operator[] (unsigned i) const
    {
        EXCEPTION("Zero-dimensional point has no data");
    }
};


#endif //_CHASTEPOINT_HPP_
