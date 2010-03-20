/*

Copyright (C) University of Oxford, 2005-2010

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
      mNodeItemWidth(0),
      mElementItemWidth(0),
      mFaceItemWidth(0),
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
      mReadContainingElementOfBoundaryElement(readContainingElementForBoundaryElements),
      mFilesAreBinary(false),
      mMeshIsHexahedral(false)
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
    std::vector<double> ret_coords(SPACE_DIM);

    // There are no attributes to read with node coordinates
    const unsigned num_attributes = 0u;
    unsigned empty = 0u;
    GetNextItemFromStream(mNodesFile, mNodesRead, ret_coords, num_attributes, empty);

    mNodesRead++;
    return ret_coords;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElementData()
{
    ElementData element_data;
    element_data.NodeIndices.resize(mNodesPerElement);
    element_data.AttributeValue = 0; // If an attribute is not read this stays as zero, otherwise overwritten.
    GetNextItemFromStream(mElementsFile, mElementsRead, element_data.NodeIndices, mNumElementAttributes,
                          element_data.AttributeValue);

    EnsureIndexingFromZero(element_data.NodeIndices);

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
        ret_indices.resize(mNodesPerBoundaryElement);

        assert(ELEMENT_DIM != 0); //Covered in earlier exception, but needed in loop guard here.
        do
        {
            face_data.AttributeValue = 1u; // If an attribute is not read this stays as one, otherwise overwritten.


            if (mReadContainingElementOfBoundaryElement)
            {
                assert(mNumFaceAttributes == 0);
                GetNextItemFromStream(mFacesFile, mFacesRead, ret_indices, 1,
                                      face_data.ContainingElement);
            }
            else
            {
                GetNextItemFromStream(mFacesFile, mFacesRead, ret_indices, mNumFaceAttributes,
                                      face_data.AttributeValue);
            }

            EnsureIndexingFromZero(ret_indices);

            mFacesRead++;
        }
        while (ELEMENT_DIM==2 && face_data.AttributeValue==0); //In triangles format we ignore internal edges (which are marked with attribute 0)
    }

    mBoundaryFacesRead++;
    face_data.NodeIndices = ret_indices;
    return face_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index)
{
    if (!mFilesAreBinary)
    {
        EXCEPTION("Random access is only implemented in mesh readers for binary mesh files.");
    }
    if (index >= mNumNodes)
    {
        EXCEPTION("Node does not exist - not enough nodes.");
    }
    // Put the file stream pointer to the right location
    mNodesFile.seekg(mNodeFileDataStart + mNodeItemWidth*index, std::ios_base::beg);
    // Read the next item.
    return GetNextNode();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetElementData(unsigned index)
{
    if (!mFilesAreBinary)
    {
        EXCEPTION("Random access is only implemented in mesh readers for binary mesh files.");
    }
    if (index >=mNumElements)
    {
        EXCEPTION("Element does not exist - not enough elements.");
    }

    // Put the file stream pointer to the right location
    mElementsFile.seekg(mElementFileDataStart + mElementItemWidth*index, std::ios_base::beg);
    // Read the next item.
    return GetNextElementData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetFaceData(unsigned index)
{
    if (!mFilesAreBinary)
    {
        EXCEPTION("Random access is only implemented in mesh readers for binary mesh files.");
    }
    if (index >=mNumFaces)
    {
        EXCEPTION("Face does not exist - not enough faces.");
    }
    // Put the file stream pointer to the right location
    mFacesFile.seekg(mFaceFileDataStart + mFaceItemWidth*index, std::ios_base::beg);
    // Read the next item.
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
    /*
     *  Reading node file header
     */
    std::string buffer;
    GetNextLineFromStream(mNodesFile, buffer);
    std::stringstream node_header_line(buffer);
    unsigned dimension;
    node_header_line >> mNumNodes >> dimension >> mNumNodeAttributes >> mMaxNodeBdyMarker;
    if (SPACE_DIM != dimension)
    {
        EXCEPTION("SPACE_DIM  != dimension read from file ");
    }
    //Is there anything else on the header line?
    std::string extras;
    node_header_line >> extras;
    if (extras == "BIN")
    {
        mFilesAreBinary = true;
        mNodeFileDataStart = mNodesFile.tellg(); // Record the position of the first byte after the header.
        mNodeItemWidth = SPACE_DIM * sizeof(double);
        //We enforce that all binary files (written by Chaste) are indexed from zero
        mIndexFromZero = true;
    }
    else
    {
        // Get the next line to see if it is indexed from zero or not
        GetNextLineFromStream(mNodesFile, buffer);
        std::stringstream node_first_line(buffer);
        unsigned first_index;
        node_first_line >> first_index;
        assert(first_index == 0 || first_index == 1);
        mIndexFromZero = (first_index == 0);
        // Close, reopen, skip header
        mNodesFile.close();
        OpenNodeFile();
        GetNextLineFromStream(mNodesFile, buffer);
    }

    /*
     *  Reading element file header
     */
    GetNextLineFromStream(mElementsFile, buffer);
    std::stringstream element_header_line(buffer);

    unsigned extra_attributes = 0;

    if (ELEMENT_DIM == SPACE_DIM)
    {
        element_header_line >> mNumElements >> mNumElementNodes >> mNumElementAttributes;

        extra_attributes = mNumElementAttributes;

        //Is there anything else on the header line?
        std::string element_extras;
        element_header_line >> element_extras;
        if (element_extras == "BIN")
        {
            //Double check for binaryness
            assert (mFilesAreBinary);
        }
        else if (element_extras == "HEX")
        {
            mMeshIsHexahedral = true;
            if ( ELEMENT_DIM == 2 )
            {
                mNodesPerElement = 4;
                mNodesPerBoundaryElement = 2;
            }
            if ( ELEMENT_DIM == 3 )
            {
                mNodesPerElement = 8;
                mNodesPerBoundaryElement = 4;
            }
        }
        else
        {
            assert (element_extras == "");
        }

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
        element_header_line >> mNumElements >> mNumFaceAttributes;

        extra_attributes = mNumFaceAttributes;

        //Is there anything else on the header line?
        std::string element_extras;
        element_header_line >> element_extras;
        if (element_extras == "BIN")
        {
            //Double check for binaryness
            assert (mFilesAreBinary);
        }

        mNodesPerElement = ELEMENT_DIM+1;
    }


    if (mFilesAreBinary)
    {
        mElementFileDataStart = mElementsFile.tellg(); // Record the position of the first byte after the header.
        mElementItemWidth = (mNodesPerElement + extra_attributes) * sizeof(unsigned);
    }

    /*
     *  Reading face/edge file header
     */
    if (ELEMENT_DIM == 1)
    {
       GetOneDimBoundary();
       mNumFaces = mOneDimBoundary.size();
    }
    else
    {
        GetNextLineFromStream(mFacesFile, buffer);
        std::stringstream face_header_line(buffer);

        face_header_line >> mNumFaces >> mNumFaceAttributes;
        assert(mNumFaceAttributes==0 || mNumFaceAttributes==1);
        // if mNumFaceAttributes=1 then loop over and set mNumFaces to be
        // the number of faces which are marked as boundary faces
        //Double check for binaryness
        std::string face_extras;
        face_header_line >> face_extras;
        assert (mFilesAreBinary == (face_extras == "BIN"));
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

    if (mFilesAreBinary)
    {
        mFaceFileDataStart = mFacesFile.tellg(); // Record the position of the first byte after the header.
        mFaceItemWidth = (ELEMENT_DIM + mNumFaceAttributes) * sizeof(unsigned);
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
template<class T>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextItemFromStream(std::ifstream& fileStream, unsigned expectedItemNumber,
                               std::vector<T>& rDataPacket, const unsigned& rNumAttributes, unsigned& rAttribute)
{
    if (mFilesAreBinary)
    {
        fileStream.read((char*)&rDataPacket[0], rDataPacket.size()*sizeof(T));
        if (rNumAttributes>0)
        {
            assert(rNumAttributes == 1);
            fileStream.read((char*) &rAttribute, sizeof(unsigned));
        }
    }
    else
    {
        std::string buffer;
        GetNextLineFromStream(fileStream,buffer);
        std::stringstream buffer_stream(buffer);

        unsigned item_index;
        buffer_stream >> item_index;

        // If we are indexing from zero our expected item number is one larger.
        expectedItemNumber += mIndexFromZero ? 0 : 1;

        if (item_index != expectedItemNumber)
        {
            std::stringstream error;
            if (!mIndexFromZero)
            { // To fix the exception message to agree with file format.
                expectedItemNumber--;
            }
            error << "Data for item " << expectedItemNumber << " missing";
            EXCEPTION(error.str());
        }

        for (unsigned i=0; i<rDataPacket.size(); i++)
        {
            buffer_stream >> rDataPacket[i];
        }

        if (rNumAttributes>0)
        {
            assert(rNumAttributes == 1);
            buffer_stream >> rAttribute;
        }
    }
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

    std::vector<unsigned> node_indices(2);
    unsigned dummy_attribute;
    //Count how many times we see each node
    std::vector<unsigned> node_count(mNumNodes);//Covers the case if it's indexed from 1
    for (unsigned element_index=0; element_index<mNumElements;element_index++)
    {
        GetNextItemFromStream(mElementsFile, element_index, node_indices, 0, dummy_attribute);
        if (!mIndexFromZero)
        {
            //Adjust so we are indexing from zero
            node_indices[0]--;
            node_indices[1]--;
        }
        node_count[node_indices[0]]++;
        node_count[node_indices[1]]++;
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
    std::string buffer;
    GetNextLineFromStream(mElementsFile, buffer);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::EnsureIndexingFromZero(std::vector<unsigned>& rNodeIndices)
{

    if (!mIndexFromZero) // If node indices do not start at zero move them all down one so they do
    {
        for (unsigned i=0; i<rNodeIndices.size(); i++)
        {
            rNodeIndices[i]--;
        }
    }

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::IsFileFormatBinary()
{
    return mFilesAreBinary;
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


/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template void TrianglesMeshReader<0,1>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, unsigned&);
template void TrianglesMeshReader<0,1>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, unsigned&);
template void TrianglesMeshReader<1,1>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, unsigned&);
template void TrianglesMeshReader<1,1>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, unsigned&);
template void TrianglesMeshReader<1,2>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, unsigned&);
template void TrianglesMeshReader<1,2>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, unsigned&);
template void TrianglesMeshReader<1,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, unsigned&);
template void TrianglesMeshReader<1,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, unsigned&);
template void TrianglesMeshReader<2,2>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, unsigned&);
template void TrianglesMeshReader<2,2>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, unsigned&);
template void TrianglesMeshReader<2,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, unsigned&);
template void TrianglesMeshReader<2,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, unsigned&);
template void TrianglesMeshReader<3,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, unsigned&);
template void TrianglesMeshReader<3,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, unsigned&);
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
