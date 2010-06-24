/*

Copyright (C) University of Oxford, 2005-2010

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

#include "AbstractChasteRegion.hpp"
#include "ChastePoint.hpp"


/**
 * This class defines a 3D cuboid and provides a method to check
 * if a given point is contained in the volume.
 */

template <unsigned SPACE_DIM>
class ChasteCuboid : public AbstractChasteRegion<SPACE_DIM>
{
private:

    /** Lower vertex of the cuboid. */
    ChastePoint<SPACE_DIM> mLowerCorner;

    /** Upper vertex of the cuboid.  The space-diagonal opposite corner of mLowerCorner. */
    ChastePoint<SPACE_DIM> mUpperCorner;

public:

    /**
     * The cuboid is defined by any of its two space-diagonal opposite corners.
     *
     * @param rLowerPoint Lower vertex of the cuboid.
     * @param rUpperPoint Upper vertex of the cuboid.
     */
    ChasteCuboid(ChastePoint<SPACE_DIM>& rLowerPoint, ChastePoint<SPACE_DIM>& rUpperPoint)
        : mLowerCorner(rLowerPoint),
          mUpperCorner(rUpperPoint)
    {
        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            if (mLowerCorner[dim] > mUpperCorner[dim])
            {
                EXCEPTION("Attempt to create a cuboid with MinCorner greater than MaxCorner in some dimension");
            }
        }
    }


    /**
     * Checks if a given 3D point is contained in the cuboid.
     *
     * @param rPointToCheck Point to be checked to be contained in the cuboid.
     */

    bool DoesContain(const ChastePoint<SPACE_DIM>& rPointToCheck) const
    {
        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            if (rPointToCheck[dim] < mLowerCorner[dim] - 100*DBL_EPSILON
                || mUpperCorner[dim] + 100* DBL_EPSILON < rPointToCheck[dim])
            {
                return false;
            }
        }
        return true;
    }

    /** @return the upper vertex of the cuboid */
    const ChastePoint<SPACE_DIM>& rGetUpperCorner() const
    {
        return mUpperCorner;
    }
    /** @return the lower vertex of the cuboid */
    const ChastePoint<SPACE_DIM>& rGetLowerCorner() const
    {
        return mLowerCorner;
    }
    
    /** 
     * @param rDimension dimension
     * @return the width in a particular dimension */
    double GetWidth(unsigned rDimension) const
    {
        assert(rDimension<SPACE_DIM);
        return mUpperCorner[rDimension] - mLowerCorner[rDimension];
    }
};

#endif /*CHASTECUBOID_HPP_*/
