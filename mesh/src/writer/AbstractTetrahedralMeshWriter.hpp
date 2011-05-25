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


#ifndef _ABSTRACTTETRAHEDRALMESHWRITER_HPP_
#define _ABSTRACTTETRAHEDRALMESHWRITER_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralMesh;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class DistributedTetrahedralMesh;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MixedDimensionMesh;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MeshWriterIterators;

#include <fstream>
#include <sstream>
#include <iostream>

#include "Exception.hpp"
#include "AbstractMeshWriter.hpp"
#include "AbstractMesh.hpp"
#include "NodeMap.hpp"





/**
 * An abstract tetrahedral mesh writer class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralMeshWriter : public AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
private:

    /**
     * Write a parallel mesh to file. Used by the serialization methods
     *
     * @param keepOriginalElementIndexing  Whether to write the mesh with the same element ordering.
     *                                     Optimisations can be applied if this is not needed. 
     */
    virtual void WriteFilesUsingParallelMesh(bool keepOriginalElementIndexing=true);

    /**
     * Create output files and add headers.
     */
    virtual void CreateFilesWithHeaders();
    
    /**
     * Append local mesh data to output files.
     */
    virtual void AppendLocalDataToFiles();
    
    /**
     * Append footers to output files.
     */
    virtual void WriteFilesFooter();

    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh; /**<Pointer to the mesh (if we are writing from the a mesh)*/

    NodeMap* mpNodeMap; /**<Node map to be used when writing a mesh that has deleted nodes*/

protected:

    unsigned mNodesPerElement; /**< Same as (ELEMENT_DIM+1), except when writing a quadratic mesh!*/
    unsigned mNodesPerBoundaryElement; /**< Same as (ELEMENT_DIM), except when writing a quadratic mesh!*/

    DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpDistributedMesh; /**< Another pointer to the mesh, produced by dynamic cast*/
    MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* mpMixedMesh; /**< Another pointer to the mesh, produced by dynamic cast*/
    MeshWriterIterators<ELEMENT_DIM,SPACE_DIM>* mpIters; /**< Handy iterators so that we know the next node/element to be written */

    bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */
    bool mWriteMetaFile; /**< Whether to write a metafile (only used by MeshylazerMeshWriter) */
    unsigned mNodeCounterForParallelMesh; /**< Used by master process for polling processes for the next node */
    unsigned mElementCounterForParallelMesh;/**< Used by master process for polling processes for the next element */
    unsigned mBoundaryElementCounterForParallelMesh;/**< Used by master process for polling processes for the next boundary element */
    unsigned mCableElementCounterForParallelMesh;/**< Used by master process for polling processes for the next cable element */
    bool mFilesAreBinary;  /**< Whether all data is to be written as binary - used in derived class TrianglesMeshWriter*/

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
     *  Destructor just deletes the node map if memory has been allocated for it
     */
    ~AbstractTetrahedralMeshWriter();

    /**
     * Write a const mesh to file. Used by the serialization methods and avoids iterators...
     *
     * @param rMesh the mesh
     * @param keepOriginalElementIndexing  Whether to write the mesh with the same element ordering.
     *                                     Optimisations can be applied if this is not needed. 
     * 
     * \todo #1322 Mesh should really be const!
     */
    virtual void WriteFilesUsingMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, 
                                     bool keepOriginalElementIndexing=true);


   /**
     * @return the coordinates of the next node to be written to file
     */
    std::vector<double> GetNextNode();


    /**
     * @return the data (indices/attributes) of the next element to be written to file
     */
    ElementData GetNextElement();

    /**
     * @return the data (indices) of the next boundary element to be written to file
     */
    ElementData GetNextBoundaryElement();

    /**
     * @return the data (indices/attributes) of the next cable element to be written to file
     */
    ElementData GetNextCableElement();
};

#endif //_ABSTRACTTETRAHEDRALMESHWRITER_HPP_
