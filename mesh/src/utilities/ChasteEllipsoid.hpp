/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef CHASTEELLIPSOID_HPP_
#define CHASTEELLIPSOID_HPP_

#include "AbstractChasteRegion.hpp"
#include "ChastePoint.hpp"


/**
 * This class defines a 3D ellipsoid and provides a method to check
 * if a given point is contained in the volume.
 */
template <unsigned SPACE_DIM>
class ChasteEllipsoid : public AbstractChasteRegion<SPACE_DIM>
{
private:
    /** Centre of the ellipsoid. */
    ChastePoint<SPACE_DIM> mCentre;

    /** Radii of the ellipsoid. */
    ChastePoint<SPACE_DIM> mRadii;

public:
    /**
     * The (axis aligned) ellipsoid is defined by its centre and its radii in the x, y and z directions.
     *
     * @param rCentre Centre of the ellipsoid.
     * @param rRadii Radii of the ellipsoid.
     */
    ChasteEllipsoid(ChastePoint<SPACE_DIM>& rCentre, ChastePoint<SPACE_DIM>& rRadii);

    /**
     * Checks if a given point is contained in the ellipsoid.
     *
     * @param rPointToCheck Point to be checked to be contained in the ellipsoid.
     */
    bool DoesContain(const ChastePoint<SPACE_DIM>& rPointToCheck) const;

    /** @return centre of the ellipsoid */
    const ChastePoint<SPACE_DIM>& rGetCentre() const;
    
    /** @return radii of the ellipsoid */
    const ChastePoint<SPACE_DIM>& rGetRadii();

};

#endif /*CHASTEELLIPSOID_HPP_*/
