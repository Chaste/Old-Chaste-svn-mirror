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


#ifndef _ABSTRACTTETRAHEDRALMESHWRITER_HPP_
#define _ABSTRACTTETRAHEDRALMESHWRITER_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralMesh;

#include <fstream>
#include <sstream>
#include <iostream>

#include "Exception.hpp"
#include "AbstractMeshWriter.hpp"
#include "NodeMap.hpp"

/**
 * An abstract tetrahedral mesh writer class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralMeshWriter : public AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
protected:

    ///\todo The following three members aren't used anywhere - remove them? (#1076)
    std::vector< std::vector<double> >::iterator mpNodeIterator; /**< Is an iterator for the node data */
    std::vector< std::vector<unsigned> >::iterator mpElementIterator; /**< Is an iterator for the element data */
    std::vector< std::vector<unsigned> >::iterator mpBoundaryFaceIterator; /**< Is an iterator for the boundary face data */

    bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */
    bool mWriteMetaFile; /**< Whether to write a metafile (only used by MeshylazerMeshWriter) */

public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param clearOutputDir  whether to clean the directory (defaults to true)
     */
    AbstractTetrahedralMeshWriter(const std::string& rDirectory,
                       const std::string& rBaseName,
                       const bool clearOutputDir=true);

    /**
     * Write a mesh to file.
     *
     * @param rMesh the mesh
     */
    void WriteFilesUsingMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     * Write a const mesh to file. Used by the serialization methods and avoids iterators...
     *
     * \todo This is a very smelly method which has copied code from the above method...
     *
     * @param rMesh the mesh
     */
    void WriteFilesUsingMesh(const AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

};

#endif //_ABSTRACTTETRAHEDRALMESHWRITER_HPP_
