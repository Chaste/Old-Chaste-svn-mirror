/*

Copyright (C) University of Oxford, 2008

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

#endif /*CHASTECUBOID_HPP_*/
