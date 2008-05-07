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


#ifndef _ABSTRACTMESHWRITER_HPP_
#define _ABSTRACTMESHWRITER_HPP_

#include "ConformingTetrahedralMesh.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractMeshReader.hpp"


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
    bool mWriteMetaFile;
public:
    /** Constructor */
    AbstractMeshWriter(const std::string &rDirectory,
                       const std::string &rBaseName,
                       const bool clearOutputDir=true)
            : mBaseName(rBaseName)
    {
        mpOutputFileHandler = new OutputFileHandler(rDirectory, clearOutputDir);
    }
    /** Destructor */
    virtual ~AbstractMeshWriter()
    {
        delete mpOutputFileHandler;
    }
    std::string GetOutputDirectory(void);
    
    void SetNextNode(std::vector<double> nextNode);
    void SetNextElement(std::vector<unsigned> nextElement);
    void SetNextBoundaryFace(std::vector<unsigned> nextFace);
    void SetNextBoundaryEdge(std::vector<unsigned> nextEdge);
    virtual void WriteFiles()=0;
    unsigned GetNumNodes()
    {
        return mNodeData.size();
    }
    unsigned GetNumElements()
    {
        return mElementData.size();
    }
    unsigned GetNumBoundaryFaces()
    {
        return mBoundaryFaceData.size();
    }
    unsigned GetNumBoundaryEdges()
    {
        return mBoundaryFaceData.size();
    }
    void WriteFilesUsingMesh(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);
    void WriteFilesUsingMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader);
};




/**
 * Return the full path to the directory where meshes will be written.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputDirectory(void)
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
     ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
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
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter =
        rMesh.GetElementIteratorBegin();
        
    while (iter != rMesh.GetElementIteratorEnd())
    {
        if ((*iter)->IsDeleted() == false)
        {
            std::vector<unsigned> indices(ELEMENT_DIM+1);
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                unsigned old_index=(*iter)->GetNodeGlobalIndex(j);
                indices[j] = node_map.GetNewIndex(old_index);
            }
            SetNextElement(indices);
        }
        
        iter++;
    }
    
    // Get a iterator over the boundary elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator boundary_iter =
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
        SetNextElement(rMeshReader.GetNextElement());
    }
    for (unsigned i=0; i<rMeshReader.GetNumFaces();i++)
    {
        SetNextBoundaryFace(rMeshReader.GetNextFace());
    }
    WriteFiles();
}

#endif //_ABSTRACTMESHWRITER_HPP_
