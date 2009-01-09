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
#ifndef VERTEXMESHWRITER2D_HPP_
#define VERTEXMESHWRITER2D_HPP_

#include "OutputFileHandler.hpp"
#include <iomanip>

/**
 * A mesh writer class for 2D vertex-based meshes.
 */
class VertexMeshWriter2d 
{
private:

    /** Output file handler member. */
    OutputFileHandler *mpOutputFileHandler;
    
    /** Base name for results files. */
    std::string mBaseName;

public:

    /**
     * Constructor.
     * 
     * @param rDirectory reference to the output directory, relative to where Chaste output is stored
     * @param rBaseName reference to the base name for results files
     * @param clearOutputDir whether to clear the output directory prior to writing files
     */
    VertexMeshWriter2d(const std::string &rDirectory,
                       const std::string &rBaseName,
                       const bool clearOutputDir=true);

    /**
     * Destructor.
     */
    virtual ~VertexMeshWriter2d();

    /**
     * Write files.
     * 
     * @param rMesh reference to the vertex-based mesh
     */
    void WriteFiles(VertexMesh<2,2>& rMesh);
};

#endif /*VERTEXMESHWRITER2D_HPP_*/
