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
#ifndef DISTANCEMAPCALCULATOR_HPP_
#define DISTANCEMAPCALCULATOR_HPP_

#include <vector>

#include "UblasIncludes.hpp"
#include "TetrahedralMesh.hpp"

/**
 * This class provides functionalities to compute a distance map in a given mesh
 * from a given surface, specifying the distance from each node to the surface.
 * 
 * The mesh is specified in the constructor, and the ComputeDistanceMap computes
 * (and returns by reference) the map.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class DistanceMapCalculator
{
private:

    /** The mesh*/
    TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& mrMesh;
    /** Number of nodes in the mesh*/
    unsigned mNumNodes;

    /**
     *  Computes the euclidean distance of two given points
     *
     *  @param pointA First point
     *  @param pointB Second point
     */
    inline double EuclideanDistanceTwoPoints(const c_vector<double, SPACE_DIM>& pointA,
                                             const c_vector<double, SPACE_DIM>& pointB) const;

    /**
     *  Given a cartesian distance, computes the associated euclidean distance
     *
     *  @param cartDistance Cartesian distance of a given point
     */
    inline double CartToEucliDistance(c_vector<double, SPACE_DIM>& cartDistance) const;


public:

    /**
     * Constructor
     * 
     * @param rMesh the mesh to compute maps for
     */
    DistanceMapCalculator(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     *  Generates a distance map of all the nodes of the mesh to the given surface
     *
     *  @param rOriginSurface set of node indexes defining the surface
     *  @param rNodeDistances distance map computed. The method will resize it if it's not big enough.
     */
    void ComputeDistanceMap(const std::vector<unsigned>& rOriginSurface,
                            std::vector<double>& rNodeDistances);
};

#endif /*DISTANCEMAPCALCULATOR_HPP_*/
