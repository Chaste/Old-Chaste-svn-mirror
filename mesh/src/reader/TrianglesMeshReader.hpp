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


#ifndef _TRIANGLESMESHREADER_HPP_
#define _TRIANGLESMESHREADER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include "AbstractMeshReader.hpp"

/**
 * Concrete version of the AbstractCachedMeshReader class.
 * Once constructed the public methods of the AbstractCachedMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TrianglesMeshReader : public AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>
{

    friend class TestTrianglesMeshReader;//for testing

private:

    bool mIndexFromZero;            /**< True if input data is numbered from zero, false otherwise */

    std::string mFilesBaseName;     /**< The base name for mesh files. */

    std::ifstream mNodesFile;       /**< The nodes file for the mesh. */
    std::ifstream mElementsFile;    /**< The elements file for the mesh. */
    std::ifstream mFacesFile;       /**< The faces (edges) file for the mesh. */
    std::ifstream mNclFile;         /**< The node connectivity list file for the mesh. */

    std::streampos mNodeFileDataStart; /**< The start of the binary data*/
    std::streamoff mNodeItemWidth;  /**< The number of bytes in a line of the node file*/
    std::streampos mElementFileDataStart; /**< The start of the binary element data*/
    std::streamoff mElementItemWidth;  /**< The number of bytes in a line of the element file*/
    std::streampos mFaceFileDataStart; /**< The start of the binary face data*/
    std::streamoff mFaceItemWidth;  /**< The number of bytes in a line of the face file*/
    std::streampos mNclFileDataStart; /**< The start of the binary data*/
    std::streamoff mNclItemWidth;  /**< The number of bytes in a line of the node file*/

    unsigned mNumNodes;             /**< Number of nodes in the mesh. */
    unsigned mNumElements;          /**< Number of elements in the mesh. */
    unsigned mNumFaces;             /**< Number of faces in the mesh. */

    unsigned mNodesRead;            /**< Number of nodes read in. */
    unsigned mElementsRead;         /**< Number of elements read in. */
    unsigned mFacesRead;            /**< Number of faces read in. */
    unsigned mBoundaryFacesRead;    /**< Number of boundary faces read in. */
    std::vector<unsigned> mOneDimBoundary; /**<Indices of nodes which are at the boundary of a 1D mesh*/

    unsigned mNumNodeAttributes;    /**< Is the number of attributes stored at each node. */
    std::vector<double> mNodeAttributes; /**<Will contain the nodal attributes at each node. Cleared and re-filled at each node*/
    unsigned mMaxNodeBdyMarker;     /**< Is the maximum node boundary marker. */
    unsigned mNumElementNodes;      /**< Is the number of nodes per element. */
    unsigned mNumElementAttributes; /**< Is the number of attributes stored for each element. */
    unsigned mNumFaceAttributes;    /**< Is the number of attributes stored for each face. */

    unsigned mOrderOfElements;      /**< The order of each element (1 for linear, 2 for quadratic). */
    unsigned mOrderOfBoundaryElements; /**< The order of each element (1 for linear, 2 for quadratic). */
    unsigned mNodesPerElement;      /**< The number of nodes contained in each element. */
    unsigned mNodesPerBoundaryElement; /**< The number of nodes in each boundary element. */

    unsigned mMaxContainingElements; /**< The maximum number of elements that any node is contained in. */

    bool mEofException; /**< Set to true when end-of-file exception is thrown (for use in a try-catch) */

    bool mReadContainingElementOfBoundaryElement; /**< Whether to read containing element info for each boundary element (obtaining by doing tetgen with the -nn flag) */
    bool mFilesAreBinary; /**< Whether to read all data as binary (determined by a magic number in the node file header)*/
    bool mMeshIsHexahedral; /**< Whether the mesh is hexahedral (determined by a magic number in the element file header) */
    bool mNclFileAvailable; /**< Whether a ncl file exists */

    char* mNodeFileReadBuffer; /**< Buffer for node file read with std::ifstream */
    char* mElementFileReadBuffer; /**< Buffer for element file read with std::ifstream */
    char* mFaceFileReadBuffer; /**< Buffer for face file read with std::ifstream */

    bool mNodePermutationDefined; /**< Wheter to consider a user-defined node permutation when reading a mesh from file.*/
    std::vector<unsigned> mPermutationVector; /**< Permutation to be considered, i-th entry of the vector contains new index for original node i.*/
    std::vector<unsigned> mInversePermutationVector; /**< Permutation inverse, stored for performance reasons.*/

//    /** The containing element for each boundary element (obtaining by doing tetgen with the -nn flag).
//     *  In a std::vector rather than the struct to save space if not read.
//     */
//    std::vector<unsigned> mContainingElementsOfBoundaryElement;
//
//    unsigned mIndexIntoContainingElementsVector; /**< Which index to use when GetNextContainingElementOfBoundaryElement() is called */

public:

    /**
     * Constructor.
     *
     * @param pathBaseName  the base name of the files from which to read the mesh data
     *    (either absolute, or relative to the current directory)
     * @param orderOfElements  the order of each element: 1 for linear, 2 for quadratic (defaults to 1)
     * @param orderOfBoundaryElements the order of each boundary element: 1 for linear, 2 for quadratic (defaults to 1. May
     *  or may not be different to orderOfElements (Note tetgen with the -o2 flag creates quadratic elements but doesn't
     *  create quadratic faces, hence the need for this third parameter)
     * @param readContainingElementsForBoundaryElements Whether to read in the containing element infomation
     *  for each boundary element (in the .face file if tetgen was run with '-nn').
     */
    TrianglesMeshReader(std::string pathBaseName,
                        unsigned orderOfElements=1,
                        unsigned orderOfBoundaryElements=1,
                        bool readContainingElementsForBoundaryElements=false);

    /**
     * Destructor
     */
     ~TrianglesMeshReader();

    /** Returns the number of elements in the mesh */
    unsigned GetNumElements() const;

    /** Returns the number of nodes in the mesh */
    unsigned GetNumNodes() const;

    /** Returns the number of faces in the mesh (synonym of GetNumEdges()) */
    unsigned GetNumFaces() const;

    /** Returns the number of attributes in the mesh */
    unsigned GetNumElementAttributes() const;

    /** Returns the number of attributes in the mesh */
    unsigned GetNumFaceAttributes() const;

    /** Resets pointers to beginning*/
    void Reset();

    /** Returns a vector of the coordinates of each node in turn */
    std::vector<double> GetNextNode();

    /** Returns a vector of the nodes of each element (and any attribute infomation, if there is any) in turn */
    ElementData GetNextElementData();

    /** Returns a vector of the nodes of each face in turn (synonym of GetNextEdgeData()) */
    ElementData GetNextFaceData();


    /**
     * @return the expected order of the element file (1=linear, 2=quadratic)
     */
    unsigned GetOrderOfElements()
    {
        return mOrderOfElements;
    }
    /**
     * @return the expected order of the boundary element file (1=linear, 2=quadratic)
     */
    unsigned GetOrderOfBoundaryElements()
    {
        return mOrderOfBoundaryElements;
    }

    /**
     * @return true if the boundary element file is linear, but contains information about neighbouring elements
     */
    bool GetReadContainingElementOfBoundaryElement()
    {
        return mReadContainingElementOfBoundaryElement;
    }

    /**
     * @return the vector of node attributes
     */
    std::vector<double>& rGetNodeAttributes();

    /**
     *  Normally throws an exception.  Only implemented for tetrahedral mesh reader of binary files.
     *
     * @param index  The global node index
     * @return a vector of the coordinates of the node
     */
    std::vector<double> GetNode(unsigned index);

    /**
     *  Normally throws an exception.  Only implemented for tetrahedral mesh reader of binary files.
     *
     * @param index  The global element index
     * @return a vector of the node indices of the element (and any attribute infomation, if there is any)
     */
    ElementData GetElementData(unsigned index);

    /**
     *  Normally throws an exception.  Only implemented for tetrahedral mesh reader of binary files.
     *
     * @param index  The global face index
     * @return a vector of the node indices of the face (and any attribute/containment infomation, if there is any)
     */
    ElementData GetFaceData(unsigned index);

    /**
     *  Normally throws an exception.  When a NCL file is available, returns a list of the elements
     *  that contain the node (only available for binary files).
     *
     * @param index  The global node index
     * @return a vector of the node indices of the face (and any attribute/containment infomation, if there is any)
     */
    std::vector<unsigned> GetContainingElementIndices(unsigned index);


    /*** Returns true if reading binary files, false if reading ascii files */
    bool IsFileFormatBinary();

    /**
     * Returns true if there is a node connectivity list (NCL) file available.
     *
     * @return whether there is a node connectivity list (NCL) file available
     */
    bool HasNclFile();

    /**
     * Sets size of std:ifstream internal read buffer. Use it for tuning I/O.
     * 
     * @param bufferSize The size of the read buffer in bytes.
     */
    void SetReadBufferSize(unsigned bufferSize);

    /**
     * Sets a node permutation to use when reading in node file.
     * 
     * @param rPermutationVector Permutation vector
     */
    void SetNodePermutation(std::vector<unsigned>& rPermutationVector);
    
private:

    /** Open mesh files. */
    void OpenFiles();

    /** Open node file. \todo Change name to OpenNodesFile for consistency with OpenElementsFile and OpenFacesFile? (#991) */
    void OpenNodeFile();

    /** Open elements file. */
    void OpenElementsFile();

    /** Open faces file. */
    void OpenFacesFile();

    /** Open node connectivity list file. */
    void OpenNclFile();

    /** Read the header from each mesh file. */
    void ReadHeaders();

    /** Close mesh files. */
    void CloseFiles();

    /**
     * Read in the next line.
     *
     * @param rFileStream  The file to read from
     * @param rRawLine  Will be filled in with the next line
     */
    void GetNextLineFromStream(std::ifstream& rFileStream, std::string& rRawLine);

    /**
     * Returns the Item details from the next line via a call to GetNextLineFromStream()
     *
     * @param rFileStream  The file to read from
     * @param expectedItemNumber  To check file syntax, what item is expected to be on the next line.
     * @param rDataPacket  Assumed to be of the right size but is allowed to contain dirty data on entry.
     * @param rNumAttributes  The number of attributes per item that we expect to read. Either #mNumFaceAttributes or #mNumElementAttributes.
     * @param rAttributes  Will be filled with the attribute values if rNumAttributes > 0, otherwise empty
     */
    template<class T>
    void GetNextItemFromStream(std::ifstream& rFileStream, unsigned expectedItemNumber,
                               std::vector<T>& rDataPacket, const unsigned& rNumAttributes, std::vector<T>& rAttributes);

    /** Get method for mFilesBaseName. */
    std::string GetMeshFileBaseName();

    /** Get method specialized to 1D meshes */
    void GetOneDimBoundary();

    /**
     * Helper method to ensure we are indexing the nodes from 0
     * (some files have them indexed from 1)
     * decides according to the property mIndexFromZero
     *
     * @param rNodeIndices  The nodes we have read in.
     */
    void EnsureIndexingFromZero(std::vector<unsigned>& rNodeIndices);
    
};

#endif //_TRIANGLESMESHREADER_HPP_
