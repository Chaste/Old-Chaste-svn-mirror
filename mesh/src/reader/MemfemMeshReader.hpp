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


#ifndef _MEMFEMMESHREADER_HPP_
#define _MEMFEMMESHREADER_HPP_

/**
 * Concrete version of the AbstractMeshReader class.
 * A MemfemMeshReader takes the base name of a set of Memfem
 * mesh files (ie. the path and name of the files without the suffices).
 * Once constructed the public methods of the AbstractMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data
 */

#include "AbstractMeshReader.hpp"
#include "Exception.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MemfemMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>
{
private:
    std::vector<std::vector<double> > TokenizeStringsToDoubles(
        std::vector<std::string> rawData);

    std::vector<std::vector<unsigned> > TokenizeStringsToInts(
        std::vector<std::string> rawData,
        unsigned dimensionOfObject,
        bool readHeader);

public:
    MemfemMeshReader(std::string pathBaseName);
    virtual ~MemfemMeshReader();
};



/**
 * The constructor takes the base name of a set of Memfem
 * mesh files (ie. the path and name of the files without the suffices)
 * and allows the data to be queried.
 * Typical use:
 *    AbstractMeshReader *spMeshReader=new MemfemMeshReader(
 *                        "pdes/tests/meshdata/Memfem_slab");
 * Also calls the superclass AbstractMeshReader's constructor
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::MemfemMeshReader(std::string pathBaseName)
        : AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>()
{

    //Open node file and store the lines as a vector of strings (minus the comments)
    std::string nodeFileName=pathBaseName+".pts";
    this->mNodeRawData=this->GetRawDataFromFile(nodeFileName);

    /* Read single line header which is the number of nodes */
    std::stringstream node_header_stream(this->mNodeRawData[0]);
    unsigned num_nodes;
    node_header_stream >> num_nodes;

    /* All Memfem data is in 3-d. */
    if (SPACE_DIM != 3  || ELEMENT_DIM != 3)
    {
        EXCEPTION("You have asked to read non-3D data. All Memfem data is in 3D.");
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

    //Open element file and store the lines as a vector of strings (minus the comments)
    std::string elementFileName=pathBaseName+".tetras";
    this->mElementRawData=this->GetRawDataFromFile(elementFileName);

    /* Read single line header which is the number of elements   */
    std::stringstream element_header_stream(this->mElementRawData[0]);
    unsigned num_elements;
    element_header_stream >> num_elements;


    // Read the rest of the element data using TokenizeStringsToInts method
    this->mElementData = TokenizeStringsToInts(this->mElementRawData,SPACE_DIM+1, true);
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

    //Open boundary face file and store the lines as a vector of strings (minus the comments)
    std::string faceFileName=pathBaseName+".tri";
    this->mFaceRawData=this->GetRawDataFromFile(faceFileName);

    /* There is no header.
     */

    // Read the face/edge data using TokenizeStringsToInts method
    this->mFaceData = TokenizeStringsToInts(this->mFaceRawData,SPACE_DIM,false);
    this->mpFaceIterator = this->mFaceData.begin();
}


/**
 * TokenizeStringsToDoubles is specific to reading node data which came from
 * a Memfem file.
 * Each string is expected to be 3 doubles (representing x,y,z)
 * Return value is a vector where each item is a vector of doubles which represents
 * position.  Indices are implicit in the vector.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<double> > MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToDoubles(
    std::vector<std::string> rawData)
{
    std::vector<std::vector<double> > tokenized_data; // Output

    //Iterate over the lines of input
    std::vector<std::string>::iterator the_iterator;
    for ( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
    {
        std::string line_of_data=*the_iterator;
        //std::cout << line_of_data << std::endl;
        std::stringstream line_stream(line_of_data);

        if (the_iterator!=rawData.begin()) //Ignore the header string
        {
            std::vector<double> current_coords;

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

    return tokenized_data;
}


/**
 * TokenizeStringsToInts is for reading element or boundary face data which came from
 * a Memfem file.
 *  Each string is expected to be:
 *  3 or 4 node indices
 *  ( 3 indices for a face, 4 for a tetrahedron)
 *  a region marker? (if it's an element)
 *  NB: Region markers are currently ignored.
 * Return value is a vector where each item is a vector of ints which represents
 * indices of nodes.
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<unsigned> > MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToInts(
    std::vector<std::string> rawData,
    unsigned dimensionOfObject,
    bool readHeader)
{
    std::vector<std::vector<unsigned> > tokenized_data;

    std::vector<std::string>::iterator the_iterator;
    for ( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
    {
        std::string line_of_data=*the_iterator;
        std::stringstream line_stream(line_of_data);


        if ( readHeader == false || the_iterator!=rawData.begin() )
        {
            std::vector<unsigned> current_indices;

            for (unsigned i = 0; i < dimensionOfObject; i++)
            {
                unsigned item_index;
                line_stream >> item_index;
                //The nodes have been indexed from one so we need to shift the indices
                item_index -= 1;
                current_indices.push_back(item_index);
            }

            tokenized_data.push_back(current_indices);
        }

    }

    return tokenized_data;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::~MemfemMeshReader()
{}
#endif //_MEMFEMMESHREADER_HPP_
