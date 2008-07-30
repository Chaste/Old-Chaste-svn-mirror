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


/**
 * Concrete version of the AbstractMeshReader class.
 * A TrianglesMeshReader takes the base name of a set of Triangles or
 * Tetgen mesh files (ie. the path and name of the files without the suffices).
 * Once constructed the public methods of the AbstractMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data
 */
#ifndef _TRIANGLESMESHREADER_H_
#define _TRIANGLESMESHREADER_H_

#include "AbstractMeshReader.hpp"
#include "Exception.hpp"
#include <cassert>

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TrianglesMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestTrianglesMeshReader;

private:
    std::vector<std::vector<double> > TokenizeStringsToDoubles(
        std::vector<std::string> rawData);

    /**
     *  Convert strings for vectors of unsigned
     *  @param rawData the raw data. Each string should look like: index value0 value1 .. valueN marker
     *  Here marker doesn't have to be present, it is ignored unless onlyMarked=true
     *  @dimensionOfObject The number of values
     *  @onlyMarked Set this to true to look at the marker and ignore any strings
     *  for which the marker is set to zero.
     */
    std::vector<std::vector<unsigned> > TokenizeStringsToInts(
        std::vector<std::string> rawData,
        unsigned dimensionOfObject,
        bool onlyMarked=false);

    void ReadFacesAsElements(std::string pathBaseName);
    void ReadEdgesAsFaces(std::string pathBaseName);

public:
    TrianglesMeshReader(std::string pathBaseName, unsigned orderOfElements=1);

    virtual ~TrianglesMeshReader();
};





/**
 * The constructor takes the base name of a set of Triangles or
 * Tetgen mesh files (ie. the path and name of the files without the suffices)
 * and allows the data to be queried.
 * Typical use:
 *    AbstractMeshReader *spMeshReader=new TrianglesMeshReader(
 *                        "pdes/tests/meshdata/disk_522_elements");
 * Also calls the superclass AbstractMeshReader's constructor
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::TrianglesMeshReader(std::string pathBaseName, unsigned orderOfElements /*=1*/)
{
    assert(orderOfElements==1 || orderOfElements==2);
    unsigned expected_num_nodes_per_elem;
    if(orderOfElements==1)
    {
        expected_num_nodes_per_elem = ELEMENT_DIM+1;
    }
    else
    {
        assert(ELEMENT_DIM==1 || ELEMENT_DIM==2); // TODO: 3D
        
        assert(SPACE_DIM==ELEMENT_DIM);
        expected_num_nodes_per_elem = (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2;
    }
    
    //Open node file and store the lines as a vector of strings (minus the comments)
    std::string nodeFileName=pathBaseName+".node";
    this->mNodeRawData=this->GetRawDataFromFile(nodeFileName);

    /* Read single line header as at:
     * http://tetgen.berlios.de/fformats.node.html
     * http://www-2.cs.cmu.edu/~quake/triangle.node.html
     */
    std::stringstream node_header_stream(this->mNodeRawData[0]);
    unsigned num_nodes;
    unsigned dimension;
    node_header_stream >> num_nodes >> dimension >> this->mNumNodeAttributes >> this->mMaxNodeBdyMarker;

    if (SPACE_DIM != dimension)
    {
        EXCEPTION("SPACE_DIM  != dimension read from file ");
    }
    // Read the rest of the node data using TokenizeStringsToDoubles method
    this->mNodeData = TokenizeStringsToDoubles(this->mNodeRawData);
    //Initialise iterator for public GetNextNode method
    this->mpNodeIterator = this->mNodeData.begin();

    //Check that the size of the data matches the information in the header
    if (num_nodes != this->mNodeData.size())
    {
        // ignored from coverage because otherwise would have to create files
        // for a bad mesh just to test this line
        #define COVERAGE_IGNORE
        EXCEPTION("Number of nodes does not match expected number declared in header");
        #undef COVERAGE_IGNORE
    }

    if (ELEMENT_DIM < SPACE_DIM)
    {
        ReadFacesAsElements(pathBaseName);
        ReadEdgesAsFaces(pathBaseName);
        return;
    }
    // the default is to read the element file

    //Open element file and store the lines as a vector of strings (minus the comments)
    std::string elementFileName=pathBaseName+".ele";
    this->mElementRawData=this->GetRawDataFromFile(elementFileName);

    /* Read single line header as at:
    * http://tetgen.berlios.de/fformats.ele.html
    * http://www-2.cs.cmu.edu/~quake/triangle.ele.html
    */
    std::stringstream element_header_stream(this->mElementRawData[0]);
    unsigned num_elements;
    element_header_stream >> num_elements >> this->mNumElementNodes >> this->mNumElementAttributes;


    if( this->mNumElementNodes != expected_num_nodes_per_elem )
    {
        std::stringstream error;
        error << "Number of nodes per elem, " << this->mNumElementNodes << ", does not match "
              << "expected number, " << expected_num_nodes_per_elem << " (which is calculated given "
              << "the order of elements chosen, " << orderOfElements << " (1=linear, 2=quadratics)";
        EXCEPTION(error.str());
    }


    assert(this->mMaxFaceBdyMarker == 0);

    // Read the rest of the element data using TokenizeStringsToInts method

//    this->mElementData = TokenizeStringsToInts(this->mElementRawData,ELEMENT_DIM+1);
this->mElementData = TokenizeStringsToInts(this->mElementRawData, expected_num_nodes_per_elem);

    this->mpElementIterator = this->mElementData.begin();

    //Check that the size of the data matches the information in the header
    if (num_elements != this->mElementData.size())
    {
        // ignored from coverage because otherwise would have to create files
        // for a bad mesh just to test this line
        #define COVERAGE_IGNORE
        EXCEPTION("Number of elements does not match expected number declared in header");
        #undef COVERAGE_IGNORE
    }

    /* Open face file and store the lines as a vector of strings (minus the comments)
     * Note that the Triangles equivalent of a face file is an edge file.
     * The two file formats are similar.  We store edges as "faces" but the superclass
     * provides a GetNextEdge method which queries this data.
     */

    std::string face_file_name;

    if (SPACE_DIM == 2)
    {
        face_file_name=pathBaseName+".edge";
    }
    else if (SPACE_DIM == 3)
    {
        face_file_name=pathBaseName+".face";
    }
    else if (SPACE_DIM == 1)
    {
        //There is no file
        //Set the mFaceData as all the nodes.
        unsigned num_faces = this->GetNumNodes();
        for (unsigned i=0; i<num_faces; i++)
        {
            std::vector<unsigned> current_item;
            current_item.push_back(i);
            this->mFaceData.push_back(current_item);
        }
        this->mpFaceIterator = this->mFaceData.begin();
        return;
    }

    this->mFaceRawData=this->GetRawDataFromFile(face_file_name);

    /* Read single line header as at:
     * http://tetgen.berlios.de/fformats.face.html
     * http://www-2.cs.cmu.edu/~quake/triangle.edge.html
     */
    std::stringstream face_header_stream(this->mFaceRawData[0]);
    unsigned num_faces;
    face_header_stream >> num_faces >> this->mMaxFaceBdyMarker;
    assert(this->mMaxFaceBdyMarker==0 || this->mMaxFaceBdyMarker==1);

    //mNumBoundaryFaces = mNumFaces; //temporary

    // Read the rest of the face/edge data using TokenizeStringsToInts method
    bool only_marked = (this->mMaxFaceBdyMarker==1);
    this->mFaceData = TokenizeStringsToInts(this->mFaceRawData,ELEMENT_DIM, only_marked);
    this->mpFaceIterator = this->mFaceData.begin();

    //Check that the size of the data matches the information in the header
    if ((!only_marked) && (num_faces != this->mFaceData.size()) )
    {
        EXCEPTION("Number of faces does not match expected number declared in header");
    }
}

/**
 * TokenizeStringsToDoubles is specific to reading node data which came from
 * a Triangles or Tetgen file.
 * http://tetgen.berlios.de/fformats.node.html
 * http://www-2.cs.cmu.edu/~quake/triangle.node.html
 *  Each string is expected to be an index number, 2 or 3 doubles (representing x,y,z)
 *  an optional list of attributes and an optional boundary marker
 *  NB: Attributes and boundary markers are currently ignored.
 * Return value is a vector where each item is a vector of doubles which represents
 * position.  Indices are implicit in the vector.
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<double> > TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToDoubles(
    std::vector<std::string> rawData)
{
    std::vector<std::vector<double> > tokenized_data; // Output

    //Iterate over the lines of input
    std::vector<std::string>::iterator the_iterator;
    for ( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
    {
        std::string line_of_data=*the_iterator;
        std::stringstream line_stream(line_of_data);

        if (the_iterator!=rawData.begin()) //Ignore the header string
        {
            std::vector<double> current_coords;
            unsigned item_number;

            //Read item index
            if (line_stream >> item_number)
            {

                //Watch for item zero (this will only happen at most once!)
                if (item_number == 0)
                {
                    this->mIndexFromZero = true;
                }

                //Throw an error if the item numbers are out of order
                if (this->mIndexFromZero == true)
                {
                    if (item_number != tokenized_data.size())
                    {
                        // ignored from coverage because otherwise would have to create files
                        // for a bad mesh just to test this line - note that it is tested for
                        // below in the case that the nodes are indexed from 1
                        #define COVERAGE_IGNORE
                        std::stringstream item_number_as_string;
                        item_number_as_string << item_number;
                        EXCEPTION("Node number " + item_number_as_string.str() + " is out of order in input file.");
                        #undef COVERAGE_IGNORE
                    }
                }
                else
                {
                    if (item_number-1 != tokenized_data.size())
                    {
                        std::stringstream item_number_as_string;
                        item_number_as_string << item_number;
                        EXCEPTION("Node number " + item_number_as_string.str() + " is out of order in input file.");
                    }
                }

                //Form the vector which represents the position of this item
                for (unsigned i = 0; i < SPACE_DIM; i++)
                {
                    double item_coord;
                    line_stream >> item_coord;
                    current_coords.push_back(item_coord);
                }

                //Put item onto main output vector
                tokenized_data.push_back(current_coords);
            }
        }
    }

    return tokenized_data;
}


/**
 * TokenizeStringsToInts is for reading element, face or edge data which came from
 * a Triangles or Tetgen file.
 * http://tetgen.berlios.de/fformats.ele.html
 * http://www-2.cs.cmu.edu/~quake/triangle.ele.html
 * http://tetgen.berlios.de/fformats.face.html
 * http://www-2.cs.cmu.edu/~quake/triangle.edge.html
 *  Each string is expected to be:
 *  an index number,
 *  2, 3 or 4 node indices
 *  ( In 2-D: 2 indices for an edge, 3 for a triangle)
 *  ( In 3-D: 3 indices for a face, 4 for a tetrahedron)
 *  an optional list of attributes (if it's an element)
 *  and an optional boundary marker (if it's an edge or face)
 *  NB: Attributes and boundary markers are currently ignored.
 * Return value is a vector where each item is a vector of ints which represents
 * indices of nodes.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<unsigned> > TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToInts(
    std::vector<std::string> rawData,
    unsigned dimensionOfObject,
    bool onlyMarked)
{
    std::vector< std::vector<unsigned> > tokenized_data;

    std::vector<std::string>::iterator the_iterator;
    for ( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
    {
        std::string line_of_data=*the_iterator;
        std::stringstream line_stream(line_of_data);

        if (the_iterator!=rawData.begin())
        {
            std::vector<unsigned> current_indices;
            unsigned item_number;

            if (line_stream >> item_number)
            {

                for (unsigned i=0; i < dimensionOfObject; i++)
                {
                    unsigned item_index;
                    line_stream >> item_index;
                    //If the nodes have been indexed from one then we need to shift the indices
                    if (this->mIndexFromZero == false)
                    {
                        item_index -= 1;
                    }
                    current_indices.push_back(item_index);
                }

                // if asked for, look at the marker and ignore if the marker
                // is equal to 0
                if(onlyMarked)
                {
                    // read the marker
                    unsigned boundary_marker;
                    line_stream >> boundary_marker;
                    assert(boundary_marker==0 || boundary_marker==1);
                    if (boundary_marker != 0)
                    {
                        tokenized_data.push_back(current_indices);
                    }
                }
                else
                {
                    tokenized_data.push_back(current_indices);
                }
            }
        }
    }

    return tokenized_data;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadFacesAsElements(std::string pathBaseName)
{
    std::string face_file_name;

    if (SPACE_DIM == 2 && ELEMENT_DIM == 1)
    {
        face_file_name=pathBaseName+".edge";
    }
    else if (SPACE_DIM == 3 && ELEMENT_DIM == 2)
    {
        face_file_name=pathBaseName+".face";
    }
    else
    {
        EXCEPTION("Can't have a zero-dimensional mesh in a one-dimensional space or a one-dimensional mesh in a three-dimensional space");
    }

    this->mElementRawData=this->GetRawDataFromFile(face_file_name);

    std::stringstream face_header_stream(this->mElementRawData[0]);
    unsigned num_elements;
    face_header_stream >> num_elements >> this->mMaxFaceBdyMarker;

    // Read the rest of the element data using TokenizeStringsToInts method
    this->mElementData = TokenizeStringsToInts(this->mElementRawData,ELEMENT_DIM+1);
    this->mpElementIterator = this->mElementData.begin();
    //Check that the size of the data matches the information in the header

    if (num_elements != this->mElementData.size())
    {
        // ignored from coverage because otherwise would have to create files
        // for a bad mesh just to test this line
        #define COVERAGE_IGNORE
        EXCEPTION("Number of elements does not match expected number declared in header");
        #undef COVERAGE_IGNORE
    }
}

/**
 * This method is specific to reading a mesh of trianges in 3D
 *
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadEdgesAsFaces(std::string pathBaseName)
{
    if (SPACE_DIM == 2 && ELEMENT_DIM == 1)
    {
        //There is no file
        //Set the mFaceData as all the nodes.
        unsigned num_faces = this->GetNumNodes();
        for (unsigned i=0; i<num_faces; i++)
        {
            std::vector<unsigned> current_item;
            current_item.push_back(i);
            this->mFaceData.push_back(current_item);
        }
        this->mpFaceIterator = this->mFaceData.begin();

        return;
    }

    //Gcov is confused by this assertion
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 3 && ELEMENT_DIM == 2);
    #undef COVERAGE_IGNORE

    std::string face_file_name;
    face_file_name=pathBaseName+".edge";

    this->mFaceRawData=this->GetRawDataFromFile(face_file_name);

    /* Read single line header as at:
     * http://tetgen.berlios.de/fformats.face.html
     * http://www-2.cs.cmu.edu/~quake/triangle.edge.html
     */
    std::stringstream face_header_stream(this->mFaceRawData[0]);
    unsigned num_faces;
    face_header_stream >> num_faces >> this->mMaxFaceBdyMarker;

    // Read the rest of the face/edge data using TokenizeStringsToInts method
    bool only_marked = (this->mMaxFaceBdyMarker==1);
    this->mFaceData = TokenizeStringsToInts(this->mFaceRawData,ELEMENT_DIM,only_marked);
    this->mpFaceIterator = this->mFaceData.begin();

    //Check that the size of the data matches the information in the header
    if (num_faces != this->mFaceData.size() && !only_marked)
    {
        // ignored from coverage because otherwise would have to create files
        // for a bad mesh just to test this line
        #define COVERAGE_IGNORE
        EXCEPTION("Number of faces does not match expected number declared in header");
        #undef COVERAGE_IGNORE
    }
}

/**
 * Destructor
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::~TrianglesMeshReader()
{
}

#endif //_TRIANGLESMESHREADER_H_
