/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "ChasteCuboid.hpp"

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

