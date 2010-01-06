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

#ifndef ABSTRACTCHASTEREGION_HPP_
#define ABSTRACTCHASTEREGION_HPP_

#include <cassert>
#include <climits> // Work around a boost bug - see #1024.

#include "ChastePoint.hpp"

/**
 * Abstract base class for Chaste regions.
 */

class AbstractChasteRegion
{

public:

    /**
     * Constructor
     */
    AbstractChasteRegion()
    {}
    
    /**
     * Checks whether the Chaste point is contained in the region. implemented in the concrete classes
     * 
     * @param rPointToCheck Point to be checked to be contained in the region
     * @return true if the point is contained, false otherwise
     */
    template <unsigned DIM>
    bool DoesContain(const ChastePoint<DIM>& rPointToCheck) const
    {
        ///\todo: This method should be pure virtual, but it's not possible to declare pure virtual templated methods. This implementation is overwritten in the concrete classes (two at the moment). Do we really need to call this method with 1D and 2D chaste points?
        NEVER_REACHED;
        return true;
    }



};
#endif /*ABSTRACTCHASTEREGION_HPP_*/
