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

#ifndef HONEYCOMBMESHGENERATOR_HPP_
#define HONEYCOMBMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "MutableMesh.hpp"

/**
 * Honeycomb mesh generator that creates a 2D honeycomb mesh (with equal distance
 * between nodes) for use in cell-centre simulations.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
class HoneycombMeshGenerator
{
protected:

    /** A pointer to the mesh this class creates */
    MutableMesh<2,2>* mpMesh;

    /** The indices of the nodes in this mesh which are 'ghost nodes'  */
    std::set<unsigned> mGhostNodeIndices;

    /** The mesh is generated by writing out a series of nodes and reading them in from this file*/
    std::string mMeshFilename;

    /** The (x) width of the domain to be constructed*/
    double mDomainWidth;

    /** The (y) depth of the domain to be constructed*/
    double mDomainDepth;

    /** The y coordinate of the bottom row of cells (ghosts if requested) */
    double mBottom;

    /** The y coordinate of the top row of cells (ghosts if requested) */
    double mTop;

    /** The number of columns of cells to put across the x coordinate of the mesh */
    unsigned mNumCellWidth;

    /** The number of rows of cells to put up the y coordinate of the mesh */
    unsigned mNumCellLength;

public:

    /**
     * Default constructor.
     *
     * @param numNodesAlongWidth  The number of cells you want alopng the bottom of the domain
     * @param numNodesAlongLength  The number of cells you want sides of the domain
     * @param ghosts  The thickness of ghost nodes to put around the edge (defaults to 0)
     * @param scaleFactor  The scale factor for the width (circumference) of the cells (defaults to 1.0)
     */
    HoneycombMeshGenerator(unsigned numNodesAlongWidth, unsigned numNodesAlongLength, unsigned ghosts=0, double scaleFactor=1.0);

    /**
     * Null constructor for derived classes to call.
     */
    HoneycombMeshGenerator()
    {
    }

    /**
     * Destructor - deletes the mesh object and pointer
     */
    virtual ~HoneycombMeshGenerator();

    /**
     * @return a 2D honeycomb mesh based on a 2D plane
     */
    virtual MutableMesh<2,2>* GetMesh();

    /**
     * Returns the indices of the nodes in the mesh which correspond to
     * real cells. This information is needed when constructing
     * a MeshBasedCellPopulationWithGhostNodes.
     *
     * @return indices of nodes
     */
    std::vector<unsigned> GetCellLocationIndices();

    /**
     * @param radius the radius of the circular mesh
     * @return a honeycomb mesh constructed to be roughly circular.
     */
    MutableMesh<2,2>* GetCircularMesh(double radius);

    /**
     * @return mDomainDepth
     */
    double GetDomainDepth();

    /**
     * @return mDomainWidth
     */
    double GetDomainWidth();
};

#endif /*HONEYCOMBMESHGENERATOR_HPP_*/
