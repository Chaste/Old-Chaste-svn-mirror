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

#ifndef CHASTECUBOID_HPP_
#define CHASTECUBOID_HPP_

#include "ChastePoint.hpp"

/**
 * This class defines a 3D cuboid and provides a method to check 
 * if a given point is contained in the volume.
 */
class ChasteCuboid
{
private:
    ChastePoint<3> mrPointA;
    ChastePoint<3> mrPointB;
    
public:
    /**
     * The cuboid is defined by any of its two space-diagonal opposite corners.
     * 
     * @param rPointA Any vertex of the cuboid.
     * @param rPointB The space-diagonal opposite corner of pointA.
     */   
    ChasteCuboid(ChastePoint<3>& rPointA, ChastePoint<3>& rPointB);
        
    /**
     * Checks if a given point is contained in the cuboid.
     * 
     * @param rPointToCheck Point to be checked to be contained in the cuboid.
     */
    bool DoesContain(const ChastePoint<3>& rPointToCheck);
};


ChasteCuboid::ChasteCuboid(ChastePoint<3>& rPointA, ChastePoint<3>& rPointB) : mrPointA(rPointA), mrPointB(rPointB)
{
    for (unsigned dim=0; dim<3; dim++)
    {
        if (mrPointA[dim] > mrPointB[dim])
        {
            EXCEPTION("Attempt to create a cuboid with MinCorner greater than MaxCorner in some dimension");
        }
    }
}

bool ChasteCuboid::DoesContain(const ChastePoint<3>& rPointToCheck)
{
    bool inside=true;
    for (unsigned dim=0; dim<3; dim++)
    {
        if (rPointToCheck[dim] < mrPointA[dim] - 100*DBL_EPSILON 
            || mrPointB[dim] + 100* DBL_EPSILON < rPointToCheck[dim])
        {
            inside=false;
            break;
        }
    }
    return inside;
}


#endif /*CHASTECUBOID_HPP_*/
