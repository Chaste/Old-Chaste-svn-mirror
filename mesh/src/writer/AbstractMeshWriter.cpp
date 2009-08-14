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

#include "AbstractMeshWriter.hpp"
#include <cassert>

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
std::string AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputDirectory()
{
    return mpOutputFileHandler->GetOutputDirectoryFullPath();
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
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextNode(std::vector<double> nextNode)
{
    assert(nextNode.size() == SPACE_DIM);
    mNodeData.push_back(nextNode);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextElement(std::vector<unsigned> nextElement)
{
    ///\todo The assertion below is only true for tetrahedral meshes - move to AbstractTetrahedralMeshWriter method?
//    assert(nextElement.size() == ELEMENT_DIM+1);
    mElementData.push_back(nextElement);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextBoundaryFace(std::vector<unsigned> nextFace)
{
    assert(nextFace.size() == ELEMENT_DIM);
    mBoundaryFaceData.push_back(nextFace);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)
{
    for (unsigned i=0; i<rMeshReader.GetNumNodes(); i++)
    {
        SetNextNode(rMeshReader.GetNextNode());
    }
    for (unsigned i=0; i<rMeshReader.GetNumElements(); i++)
    {
        SetNextElement(rMeshReader.GetNextElementData().NodeIndices);
    }
    for (unsigned i=0; i<rMeshReader.GetNumFaces(); i++)
    {
        SetNextBoundaryFace(rMeshReader.GetNextFaceData().NodeIndices);
    }
    WriteFiles();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
    const std::vector<unsigned>& rNodePermutation)
{
    if (rNodePermutation.size() == 0)
    {
       WriteFilesUsingMeshReader(rMeshReader);
    }
    else
    {
        mNodeData.resize(rMeshReader.GetNumNodes());
        for (unsigned i=0; i<rMeshReader.GetNumNodes(); i++)
        {
            assert(rNodePermutation[i] < rMeshReader.GetNumNodes());
            mNodeData[ rNodePermutation[i] ] = rMeshReader.GetNextNode();
        }

        for (unsigned i=0; i<rMeshReader.GetNumElements(); i++)
        {
            ElementData element = rMeshReader.GetNextElementData();

            for (unsigned local_index=0; local_index<element.NodeIndices.size(); local_index++)
            {
                unsigned old_index = element.NodeIndices[local_index];
                element.NodeIndices[local_index] = rNodePermutation[old_index];
            }

            SetNextElement(element.NodeIndices);
        }

        for (unsigned i=0; i<rMeshReader.GetNumFaces(); i++)
        {
            ElementData face = rMeshReader.GetNextFaceData();

            for (unsigned local_index=0; local_index<face.NodeIndices.size(); local_index++)
            {
                unsigned old_index = face.NodeIndices[local_index];
                face.NodeIndices[local_index] = rNodePermutation[old_index];
            }

            SetNextBoundaryFace(face.NodeIndices);
        }
        WriteFiles();
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractMeshWriter<1,1>;
template class AbstractMeshWriter<1,2>;
template class AbstractMeshWriter<1,3>;
template class AbstractMeshWriter<2,2>;
template class AbstractMeshWriter<2,3>;
template class AbstractMeshWriter<3,3>;
