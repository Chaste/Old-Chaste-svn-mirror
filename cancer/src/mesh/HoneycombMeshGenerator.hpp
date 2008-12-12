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
#ifndef HONEYCOMBMESHGENERATOR_HPP_
#define HONEYCOMBMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "Cylindrical2dMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"

#ifndef SPECIAL_SERIAL
#include "PetscTools.hpp"
#endif //SPECIAL_SERIAL


/**
 *  Generator of honeycomb mesh
 *
 *  This class takes in options such as width, height, number of ghost nodes
 *  and generates a honeycomb (with distance between nodes=1) mesh, and ghost
 *  node info. NOTE: the user should delete the mesh after use.
 */
/// \todo Lots of these methods and members need documentation (see #736)
class HoneycombMeshGenerator
{
private:

    MutableMesh<2,2>* mpMesh;
    std::set<unsigned> mGhostNodeIndices;
    std::string mMeshFilename;
    double mCryptWidth;
    double mCryptDepth;
    double mBottom;
    double mTop;
    unsigned mNumCellWidth;
    unsigned mNumCellLength;
    bool mCylindrical;

    /**
     * Periodic honeycomb mesh maker
     */
    void Make2dPeriodicCryptMesh(double width, unsigned ghosts);

public:

    /**
     * Crypt Periodic Honeycomb Mesh Generator
     *
     * Overwritten constructor for a mesh so mesh can be compressed by changing crypt width
     *
     * @param numNodesAlongWidth  The number of stem cells you want
     * @param numNodesAlongLength  The number of cells you want along crypt axis
     * @param ghosts The thickness of ghost nodes to put around the edge (defaults to 3)
     * @param cylindrical whether the mesh should be cylindrically periodic (defaults to true)
     * @param scaleFactor  The scale factor for the width (circumference) of the cells
     *
     * Note: this class creates a cancer params instance and sets the crypt width and length
     * accordingly in the parameters class.
     */
    HoneycombMeshGenerator(unsigned numNodesAlongWidth, unsigned numNodesAlongLength, unsigned ghosts=3, bool cylindrical=true, double scaleFactor=1.0);

    /**
     * Destructor - deletes the mesh object and pointer
     */
    ~HoneycombMeshGenerator();
    
    MutableMesh<2,2>* GetMesh();

    Cylindrical2dMesh* GetCylindricalMesh();

    std::set<unsigned> GetGhostNodeIndices();

    MutableMesh<2,2>* GetCircularMesh(double radius);

};
#endif /*HONEYCOMBMESHGENERATOR_HPP_*/
