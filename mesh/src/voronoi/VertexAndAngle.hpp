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

#ifndef VERTEXANDANGLE_HPP_
#define VERTEXANDANGLE_HPP_

#include "UblasCustomFunctions.hpp"
#include <cmath>

/**
 * Helper class containing a pointer to a vertex of the face
 * and the polar angle from the centre of the face to this
 * vertex.
 *
 * \todo does this really need to be templated?
 */
template<unsigned DIM>
class VertexAndAngle
{
private:

    /** Pointer to a vertex. */
    c_vector<double, DIM>* mpVertex;

    /** Polar angle. */
    double mAngle;

public:

    /**
     * Compute the polar angle of the point (x,y)  in the interval (-PI, PI]
     * and use this to set mAngle.
     *
     * @param x x-coordinate
     * @param y y-coordinate
     */
    void ComputeAndSetAngle(double x, double y);

    /**
     * Set mpVertex.
     *
     * @param pVertex location of the vertex
     */
    void SetVertex(c_vector<double, DIM>* pVertex);

    /**
     * Less-than angle comparison operator.
     *
     * @param rOther the VertexAndAngle object to compare to
     */
    bool operator<(const VertexAndAngle& rOther) const;

    /**
     * Get method for mAngle.
     */
    double GetAngle();

    /**
     * Get method for mpVertex.
     */
    c_vector<double, DIM>* GetVertex();
};

#endif /*VERTEXANDANGLE_HPP_*/
