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


#ifndef _ABSTRACTMESHWRITER_HPP_
#define _ABSTRACTMESHWRITER_HPP_

#include "AbstractMesh.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractMeshReader.hpp"
#include "NodeMap.hpp"

/**
 * An abstract mesh writer class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractMeshWriter
{
protected:

    OutputFileHandler *mpOutputFileHandler; /**< Output file handler */
    std::string mBaseName; /**< Base name for the input files */

    std::vector< std::vector<double> > mNodeData; /**< Is an array of node coordinates ((i,j)th entry is the jth coordinate of node i)*/
    std::vector< std::vector<unsigned> > mElementData; /**< Is an array of the nodes in each element ((i,j)th entry is the jth node of element i) */
    std::vector< std::vector<unsigned> > mBoundaryFaceData; /**< Is an array of the nodes on each boundary face ((i,j)th entry is the jth node of face i) */

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
    AbstractMeshWriter(const std::string &rDirectory,
                       const std::string &rBaseName,
                       const bool clearOutputDir=true);

    /**
     * Destructor.
     */
    virtual ~AbstractMeshWriter();

    /**
     * Return the full path to the directory where meshes will be written.
     */
    std::string GetOutputDirectory();

    /**
     * Add an entry to mNodeData.
     * 
     * @param nextNode coordinates of the node to add
     */
    void SetNextNode(std::vector<double> nextNode);

    /**
     * Add an entry to mElementData.
     * 
     * @param nextElement array of the nodes in the element to add
     */
    virtual void SetNextElement(std::vector<unsigned> nextElement);

    /**
     * Add an entry to mBoundaryFaceData.
     * 
     * @param nextFace array of the nodes on the boundary face to add
     */
    void SetNextBoundaryFace(std::vector<unsigned> nextFace);

    /**
     * Add an entry to mBoundaryFaceData.
     * 
     * @param nextEdge array of the nodes on the boundary face to add
     */
    void SetNextBoundaryEdge(std::vector<unsigned> nextEdge);

    /**
     * Write mesh data to files.
     * This method must be overridden in concrete classes.
     */
    virtual void WriteFiles()=0;

    /**
     * Get the number of nodes in the mesh.
     */
    unsigned GetNumNodes();

    /**
     * Get the number of elements in the mesh.
     */
    unsigned GetNumElements();

    /**
     * Get the number of boundary elements in the mesh.
     */
    unsigned GetNumBoundaryFaces();

    /**
     * Get the number of boundary faces in the mesh.
     */
    unsigned GetNumBoundaryEdges();

    /**
     * Write a mesh to file.
     * 
     * @param rMesh the mesh
     */
    void WriteFilesUsingMesh(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     * Read in a mesh and write it to file.
     * 
     * @param rMeshReader the mesh reader
     */
    void WriteFilesUsingMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader);

    /**
     * Read in a mesh and a given permutation of the node indices, and write the permuted mesh to file.
     * 
     * @param rMeshReader the mesh reader
     * @param rNodePermutation the node permutation
     */
    void WriteFilesUsingMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                   std::vector<unsigned>& rNodePermutation);
};


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::AbstractMeshWriter(const std::string &rDirectory,
                   const std::string &rBaseName,
                   const bool clearOutputDir)
    : mBaseName(rBaseName)
{
    mpOutputFileHandler = new OutputFileHandler(rDirectory, clearOutputDir);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::~AbstractMeshWriter()
{
    delete mpOutputFileHandler;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNodeData.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mElementData.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryFaces()
{
    return mBoundaryFaceData.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryEdges()
{
    return mBoundaryFaceData.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputDirectory()
{
    return mpOutputFileHandler->GetOutputDirectoryFullPath();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextNode(std::vector<double> nextNode)
{
    assert (nextNode.size() == SPACE_DIM);
    mNodeData.push_back(nextNode);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextElement(std::vector<unsigned> nextElement)
{
    assert (nextElement.size() == ELEMENT_DIM+1);
    mElementData.push_back(nextElement);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextBoundaryFace(std::vector<unsigned> nextFace)
{
    assert (nextFace.size() == ELEMENT_DIM);
    mBoundaryFaceData.push_back(nextFace);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextBoundaryEdge(std::vector<unsigned> nextEdge)
{
    SetNextBoundaryFace(nextEdge);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(
     AbstractMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    NodeMap node_map(rMesh.GetNumAllNodes());
    unsigned new_index=0;
    for (unsigned i=0; i<(unsigned)rMesh.GetNumAllNodes();i++)
    {
        Node<SPACE_DIM>* p_node = rMesh.GetNode(i);

        if (p_node->IsDeleted() == false)
        {
            std::vector<double> coords(SPACE_DIM);
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                coords[j] = p_node->GetPoint()[j];
            }
            SetNextNode(coords);
            node_map.SetNewIndex(i,new_index++);
        }
        else
        {
            node_map.SetDeleted(i);
        }
    }
    assert(new_index==(unsigned)rMesh.GetNumNodes());

    // Get an iterator over the elements of the mesh
    typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter =
        rMesh.GetElementIteratorBegin();

    while (iter != rMesh.GetElementIteratorEnd())
    {
        if ((*iter)->IsDeleted() == false)
        {
            std::vector<unsigned> indices((*iter)->GetNumNodes());
            for (unsigned j=0; j<indices.size(); j++)
            {
                unsigned old_index=(*iter)->GetNodeGlobalIndex(j);
                indices[j] = node_map.GetNewIndex(old_index);
            }
            this->SetNextElement(indices);
        }

        iter++;
    }

    // Get a iterator over the boundary elements of the mesh
    typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator boundary_iter =
        rMesh.GetBoundaryElementIteratorBegin();
    while (boundary_iter != rMesh.GetBoundaryElementIteratorEnd())
    {
        if ((*boundary_iter)->IsDeleted() == false)
        {
            std::vector<unsigned> indices(ELEMENT_DIM);
            for (unsigned j=0; j<ELEMENT_DIM; j++)
            {
                unsigned old_index=(*boundary_iter)->GetNodeGlobalIndex(j);
                indices[j] = node_map.GetNewIndex(old_index);
            }
            SetNextBoundaryFace(indices);
        }
        boundary_iter++;
    }
    WriteFiles();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)
{
    for (unsigned i=0; i<rMeshReader.GetNumNodes();i++)
    {
        SetNextNode(rMeshReader.GetNextNode());
    }
    for (unsigned i=0; i<rMeshReader.GetNumElements();i++)
    {
        this->SetNextElement(rMeshReader.GetNextElementData().NodeIndices);
    }
    for (unsigned i=0; i<rMeshReader.GetNumFaces();i++)
    {
        this->SetNextBoundaryFace(rMeshReader.GetNextFaceData().NodeIndices);
    }
    WriteFiles();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
    std::vector<unsigned>& rNodePermutation)
{
    if (rNodePermutation.size() == 0)
    {
       WriteFilesUsingMeshReader(rMeshReader);
    }
    else
    {
        mNodeData.resize(rMeshReader.GetNumNodes());
        for (unsigned i=0; i<rMeshReader.GetNumNodes();i++)
        {
            assert(rNodePermutation[i] < rMeshReader.GetNumNodes());
            mNodeData[ rNodePermutation[i] ] = rMeshReader.GetNextNode();
        }

        for (unsigned i=0; i<rMeshReader.GetNumElements();i++)
        {
            ElementData element = rMeshReader.GetNextElementData();

            for(unsigned local_index=0; local_index<element.NodeIndices.size(); local_index++)
            {
                unsigned old_index = element.NodeIndices[local_index];
                element.NodeIndices[local_index] = rNodePermutation[old_index];
            }

            SetNextElement(element.NodeIndices);
        }

        for (unsigned i=0; i<rMeshReader.GetNumFaces();i++)
        {
            ElementData face = rMeshReader.GetNextFaceData();

            for(unsigned local_index=0; local_index<face.NodeIndices.size(); local_index++)
            {
                unsigned old_index = face.NodeIndices[local_index];
                face.NodeIndices[local_index] = rNodePermutation[old_index];
            }

            SetNextBoundaryFace(face.NodeIndices);
        }
        WriteFiles();
    }
}

#endif //_ABSTRACTMESHWRITER_HPP_
