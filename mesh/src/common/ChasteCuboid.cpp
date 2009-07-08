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

#include "ChasteCuboid.hpp"
#include "Exception.hpp"

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

bool ChasteCuboid::DoesContain(const ChastePoint<3U>& rPointToCheck)
{
    bool inside = true;
    for (unsigned dim=0; dim<3; dim++)
    {
        if (rPointToCheck[dim] < mrPointA[dim] - 100*DBL_EPSILON
            || mrPointB[dim] + 100* DBL_EPSILON < rPointToCheck[dim])
        {
            inside = false;
            break;
        }
    }
    return inside;
}

bool ChasteCuboid::DoesContain(const ChastePoint<2U>& rPointToCheck)
{
//    EXCEPTION("Wrong argument type in ChasteCuboid::DoesContain. Can only be used in 3D");
    bool inside = true;
    for (unsigned dim=0; dim<2; dim++)
    {
        if (rPointToCheck[dim] < mrPointA[dim] - 100*DBL_EPSILON
            || mrPointB[dim] + 100* DBL_EPSILON < rPointToCheck[dim])
        {
            inside = false;
            break;
        }
    }
    return inside;
}

bool ChasteCuboid::DoesContain(const ChastePoint<1U>& rPointToCheck)
{
//    EXCEPTION("Wrong argument type in ChasteCuboid::DoesContain. Can only be used in 3D");
    bool inside = true;
    for (unsigned dim=0; dim<1; dim++)
    {
        if (rPointToCheck[dim] < mrPointA[dim] - 100*DBL_EPSILON
            || mrPointB[dim] + 100* DBL_EPSILON < rPointToCheck[dim])
        {
            inside = false;
            break;
        }
    }
    return inside;
}
