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

#include "TrianglesMeshReader.hpp"
#include "Exception.hpp"
#include <cassert>
#include <sstream>
#include <iostream>

const static char* NODES_FILE_EXTENSION = ".node";
const static char* ELEMENTS_FILE_EXTENSION = ".ele";
const static char* FACES_FILE_EXTENSION = ".face";
const static char* EDGES_FILE_EXTENSION = ".edge";

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::TrianglesMeshReader(std::string pathBaseName, 
                                                                 unsigned orderOfElements, 
                                                                 unsigned orderOfBoundaryElements,
                                                                 bool readContainingElementForBoundaryElements)
    : mFilesBaseName(pathBaseName),
      mNumNodes(0),
      mNumElements(0),
      mNumFaces(0),
      mNodesRead(0),
      mElementsRead(0),
      mFacesRead(0),
      mBoundaryFacesRead(0),
      mNumElementAttributes(0),
      mNumFaceAttributes(0),
      mOrderOfElements(orderOfElements),
      mOrderOfBoundaryElements(orderOfBoundaryElements),
      mEofException(false),
      mReadContainingElementOfBoundaryElement(readContainingElementForBoundaryElements)
{
    // Only linear and quadratic elements
    assert(orderOfElements==1 || orderOfElements==2);
    if ( mOrderOfBoundaryElements == 2 &&  mReadContainingElementOfBoundaryElement)
    {
        EXCEPTION("Boundary element file should not have containing element info if it is quadratic");
    }
    if (mOrderOfElements==1)
    {
        mNodesPerElement = ELEMENT_DIM+1;
    }
    else
    {
        #define COVERAGE_IGNORE
        assert(SPACE_DIM==ELEMENT_DIM);
        #undef COVERAGE_IGNORE
        mNodesPerElement = (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2;
    }

    if (mOrderOfBoundaryElements==1)
    {
        mNodesPerBoundaryElement = ELEMENT_DIM;
    }
    else
    {
        #define COVERAGE_IGNORE
        assert(SPACE_DIM==ELEMENT_DIM);
        #undef COVERAGE_IGNORE
        mNodesPerBoundaryElement = ELEMENT_DIM*(ELEMENT_DIM+1)/2;
    }

    mIndexFromZero = false; // Initially assume that nodes are not numbered from zero

    OpenFiles();
    ReadHeaders();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mNumFaces;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumEdges() const
{
    return mNumFaces;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElementAttributes() const
{
    return mNumElementAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaceAttributes() const
{
    return mNumFaceAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::Reset()
{
    CloseFiles();

    mNodesRead = 0;
    mElementsRead = 0;
    mFacesRead = 0;
    mBoundaryFacesRead = 0;
    mEofException = false;

    OpenFiles();
    ReadHeaders();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    std::vector<double> ret_coords;

    std::string buffer;
    GetNextLineFromStream(mNodesFile, buffer);

    std::stringstream buffer_stream(buffer);

    unsigned index;
    buffer_stream >> index;

    unsigned offset = mIndexFromZero ? 0 : 1;
    if (index != mNodesRead+offset)
    {
        std::stringstream error;
        error << "Data for node " << mNodesRead << " missing";
        EXCEPTION(error.str());
    }

    double coord;

    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        buffer_stream >> coord;
        ret_coords.push_back(coord);
    }

    mNodesRead++;

    return ret_coords;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElementData()
{
    ElementData element_data;

    std::string buffer;
    GetNextLineFromStream(mElementsFile, buffer);

    std::stringstream buffer_stream(buffer);

    unsigned element_index;
    buffer_stream >> element_index;

    unsigned offset = mIndexFromZero ? 0 : 1;
    if (element_index != mElementsRead+offset)
    {
        std::stringstream error;
        error << "Data for element " << mElementsRead << " missing";
        EXCEPTION(error.str());
    }

    unsigned node_index;
    for (unsigned i=0; i<mNodesPerElement; i++)
    {
        buffer_stream >> node_index;
        element_data.NodeIndices.push_back(node_index - offset);
    }

    if (mNumElementAttributes > 0)
    {
        assert(mNumElementAttributes==1);

        unsigned attribute_value;
        buffer_stream >> attribute_value;
        element_data.AttributeValue = attribute_value;
    }
    else
    {
        element_data.AttributeValue = 0;
    }

    mElementsRead++;
    return element_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextFaceData()
{    
    ElementData face_data;
    std::vector<unsigned> ret_indices;

    // In the first case there's no file, all the nodes are set as faces
    if (ELEMENT_DIM == 1)
    {
        ret_indices.push_back( mOneDimBoundary[mBoundaryFacesRead] );
    }
    else
    {
        unsigned offset = mIndexFromZero ? 0 : 1;

        assert(ELEMENT_DIM != 0); //Covered in earlier exception, but needed in loop guard here.
        do
        {
            ret_indices.clear();

            std::string buffer;
            GetNextLineFromStream(mFacesFile, buffer);

            std::stringstream buffer_stream(buffer);

            unsigned face_index;
            buffer_stream >> face_index;

            if (face_index != mFacesRead+offset)
            {
                std::stringstream error;
                error << "Data for face " << mFacesRead << " missing";
                EXCEPTION(error.str());
            }

            unsigned node_index;
            for (unsigned i=0; i<mNodesPerBoundaryElement; i++)
            {
                buffer_stream >> node_index;
                ret_indices.push_back(node_index-offset);
            }

            if (mNumFaceAttributes>0)
            {
                assert(mNumFaceAttributes==1);

                unsigned attribute_value;
                buffer_stream >> attribute_value;
                face_data.AttributeValue = attribute_value;
            }
            else
            {
                face_data.AttributeValue = 1u; //Do not ignore
            }
            
            if (mReadContainingElementOfBoundaryElement)
            {
                unsigned containing_element_index;
                buffer_stream >> containing_element_index;
                
                face_data.ContainingElement = containing_element_index;
            }

            mFacesRead++;
        }
        while (ELEMENT_DIM==2 && face_data.AttributeValue==0); //In triangles format we ignore internal edges (which are marked with attribute 0)
    }

    mBoundaryFacesRead++;
    face_data.NodeIndices = ret_indices;
    return face_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextEdgeData()
{
    return GetNextFaceData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenFiles()
{
    OpenNodeFile();
    OpenElementsFile();
    OpenFacesFile();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenNodeFile()
{
    // Nodes definition
    std::string file_name = mFilesBaseName + NODES_FILE_EXTENSION;
    mNodesFile.open(file_name.c_str());
    if (!mNodesFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenElementsFile()
{
    // Elements definition
    std::string file_name;
    if (ELEMENT_DIM == SPACE_DIM)
    {
        file_name = mFilesBaseName + ELEMENTS_FILE_EXTENSION;
    }
    else
    {
        if (ELEMENT_DIM == 1)
        {
            file_name = mFilesBaseName + EDGES_FILE_EXTENSION;
        }
        else if (ELEMENT_DIM == 2)
        {
            file_name = mFilesBaseName + FACES_FILE_EXTENSION;
        }
        else
        {
            EXCEPTION("Can't have a zero-dimensional mesh in a one-dimensional space");
        }
    }

    mElementsFile.open(file_name.c_str());
    if (!mElementsFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenFacesFile()
{  
    // Faces/edges definition
    std::string file_name;
    if (ELEMENT_DIM == 3)
    {
        file_name = mFilesBaseName + FACES_FILE_EXTENSION;
    }
    else if (ELEMENT_DIM == 2)
    {
        file_name = mFilesBaseName + EDGES_FILE_EXTENSION;
    }
    else //if (ELEMENT_DIM == 1)
    {
        // There is no file, data will be read from the node file (with boundaries marked)
        return;
    }

    mFacesFile.open(file_name.c_str());
    if (!mFacesFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadHeaders()
{
    std::string buffer;
    GetNextLineFromStream(mNodesFile, buffer);
    std::stringstream buffer_stream(buffer);
    unsigned dimension;
    buffer_stream >> mNumNodes >> dimension >> mNumNodeAttributes >> mMaxNodeBdyMarker;
    if (SPACE_DIM != dimension)
    {
        EXCEPTION("SPACE_DIM  != dimension read from file ");
    }

    // Get the next line to see if it is indexed from zero or not
    GetNextLineFromStream(mNodesFile, buffer);
    std::stringstream buffer_stream_ii(buffer);

    unsigned first_index;
    buffer_stream_ii >> first_index;
    assert(first_index == 0 || first_index == 1);
    mIndexFromZero = (first_index == 0);

    // Close, reopen, skip header
    mNodesFile.close();
    OpenNodeFile();
    GetNextLineFromStream(mNodesFile, buffer);

    /// \todo: rename std::stringstream variables
    GetNextLineFromStream(mElementsFile, buffer);
    std::stringstream buffer_stream2(buffer);

    if (ELEMENT_DIM == SPACE_DIM)
    {
        buffer_stream2 >> mNumElements >> mNumElementNodes >> mNumElementAttributes;

        if ( mNumElementNodes != mNodesPerElement )
        {
            std::stringstream error;
            error << "Number of nodes per elem, " << mNumElementNodes << ", does not match "
                  << "expected number, " << mNodesPerElement << " (which is calculated given "
                  << "the order of elements chosen, " << mOrderOfElements << " (1=linear, 2=quadratics)";
            EXCEPTION(error.str());
        }
    }
    else
    {
        buffer_stream2 >> mNumElements >> mNumFaceAttributes;

        mNodesPerElement = ELEMENT_DIM+1;
    }

    if (ELEMENT_DIM == 1)
    {
       GetOneDimBoundary();
       mNumFaces = mOneDimBoundary.size();
    }
    else
    {
        GetNextLineFromStream(mFacesFile, buffer);
        std::stringstream buffer_stream3(buffer);

        buffer_stream3 >> mNumFaces >> mNumFaceAttributes;
        assert(mNumFaceAttributes==0 || mNumFaceAttributes==1);
        // if mNumFaceAttributes=1 then loop over and set mNumFaces to be
        // the number of faces which are marked as boundary faces
        if ((mNumFaceAttributes==1))
        {
            unsigned num_boundary_faces = 0;
            bool end_of_file=false;
            while (!end_of_file)
            {
                try
                {
                    GetNextFaceData();
                    num_boundary_faces++;
                }
                catch(Exception& e)
                {
                    if(mEofException)
                    {
                        end_of_file = true;
                    }
                    else
                    {
                        throw e;
                    }
                }
            }
            mNumFaces = num_boundary_faces;

//// This exception would be helpful to have until #1116 is done, unfortunately some meshes do
//// actually have no boundary elements (eg closed 2d meshes in 3d space).
//            if(mNumFaces==0)
//            {
//                EXCEPTION("No boundary elements found. NOTE: elements in face/edge file with an attribute value of 0 are considered to be internal (non-boundary) elements");
//            }

            // close the file, reopen, and skip the header again
            mFacesFile.close();
            mFacesFile.clear(); // Older versions of gcc don't explicitly reset "fail" and "eof" flags in std::ifstream after calling close()
            OpenFacesFile();
            GetNextLineFromStream(mFacesFile, buffer);
            mFacesRead = 0;
            mBoundaryFacesRead = 0;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::CloseFiles()
{
    mNodesFile.close();
    mElementsFile.close();
    mFacesFile.close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextLineFromStream(std::ifstream& fileStream, std::string& rRawLine)
{
    bool line_is_blank;
    mEofException = false;
    do
    {
        getline(fileStream, rRawLine);
        if (fileStream.eof())
        {
            mEofException = true;
            /// \todo: improve this error message
            EXCEPTION("File contains incomplete data");
        }

        // Get rid of any comment
        rRawLine = rRawLine.substr(0, rRawLine.find('#',0));

        line_is_blank = (rRawLine.find_first_not_of(" \t",0) == std::string::npos);
    }
    while (line_is_blank);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetMeshFileBaseName()
{
    return mFilesBaseName;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetOneDimBoundary()
{
    assert(ELEMENT_DIM == 1);
    if (!mOneDimBoundary.empty())
    {
        //We have already read this and have reset the reader (probably from the mesh class)...
        return;
    }
    
    std::string buffer;
    unsigned dummy, node1, node2;
    
    //Count how many times we see each node
    std::vector<unsigned> node_count(mNumNodes);//Covers the case if it's indexed from 1
    for (unsigned element_index=0; element_index<mNumElements;element_index++)
    {
        GetNextLineFromStream(mElementsFile, buffer);//Header
        std::stringstream buffer_stream(buffer);
        buffer_stream >> dummy >> node1 >> node2; //Maybe a region marker too
        if (!mIndexFromZero)
        {
            //Adjust so we are indexing from zero
            node1--;
            node2--;
        }
        node_count[node1]++;
        node_count[node2]++;
    }
    //Find the ones which are terminals (only one mention)
    for (unsigned node_index=0; node_index<mNumNodes;node_index++)
    {
        if (node_count[node_index] == 1u)
        {
            mOneDimBoundary.push_back(node_index);
        }
    }
    
    // close the file, reopen, and skip the header again
    mElementsFile.close();
    mElementsFile.clear(); // Older versions of gcc don't explicitly reset "fail" and "eof" flags in std::ifstream after calling close()
    OpenElementsFile();
    GetNextLineFromStream(mElementsFile, buffer);
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class TrianglesMeshReader<0,1>;
template class TrianglesMeshReader<1,1>;
template class TrianglesMeshReader<1,2>;
template class TrianglesMeshReader<1,3>;
template class TrianglesMeshReader<2,2>;
template class TrianglesMeshReader<2,3>;
template class TrianglesMeshReader<3,3>;
