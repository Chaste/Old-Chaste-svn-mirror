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


#ifndef VORONOICELL_HPP_
#define VORONOICELL_HPP_

#include "UblasCustomFunctions.hpp"
#include "Face.hpp"
#include <cmath>


/**
 * A VoronoiCell class for use in the VoronoiTessellation class.
 */
class VoronoiCell
{
private:

    /**
     * Return whether two faces are equal.
     * 
     * @param face1 the first face
     * @param orientation1 whether the first face is oriented
     * @param face2 the second face
     * @param orientation2 whether the second face is oriented
     */
    bool EqualFaces(Face<3>& face1, bool orientation1, Face<3>& face2, bool orientation2);

public:

    /**
     * Faces of the VoronoiCell, which should be distinct.
     * 
     * This member variable is public as it is accessed directly by 
     * VoronoiTessellation methods.
     */
    std::vector< Face<3>* > mFaces;

    /**
     * How each face is oriented. From the perspective of the centre 
     * of the VoronoiCell, the vertices of each face should be ordered 
     * clockwise. If and only if this is false, the order of vertices 
     * in the corresponding face should be reversed.
     *
     * N.B. Most faces belong to two VoronoiCell, but with opposite 
     * orientations. This allows us to reuse the face data across the 
     * two cells.
     * 
     * This member variable is public as it is accessed directly by 
     * VoronoiTessellation methods.
     */
    std::vector<bool> mOrientations;

    /**
     * The centre of the VoronoiCell. 
     * 
     * This member variable is public as it is accessed directly by 
     * VoronoiTessellation methods.
     */
    c_vector<double, 3> mCellCentre;

    /**
     * Test whether two VoronoiCells are equal.
     *
     * Two VoronoiCells are equal if their set of faces are equal 
     * (including whether the faces have the same orientations).
     * 
     * @param rOtherCell the VoronoiCell to compare to
     */
    bool operator==(VoronoiCell& rOtherCell);

    /**
     * Get the centre of the VoronoiCell.
     */
    c_vector<double, 3>& rGetVoronoiCellCentre();

};

#endif /*VORONOICELL_HPP_*/
