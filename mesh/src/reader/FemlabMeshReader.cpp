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

#include "FemlabMeshReader.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::FemlabMeshReader(std::string pathBaseName,
                                                           std::string nodeFileName,
                                                           std::string elementFileName,
                                                           std::string edgeFileName)
{

    // Open node file and store the lines as a vector of strings (minus the comments)
    nodeFileName = pathBaseName + nodeFileName;
    this->mNodeRawData = this->GetRawDataFromFile(nodeFileName);

    // Read the node data using TokenizeStringsToDoubles method
    this->mNodeData = TokenizeStringsToDoubles(this->mNodeRawData);

    // Initialise iterator for public GetNextNode method
    this->mpNodeIterator = this->mNodeData.begin();


    // Open element file and store the lines as a vector of strings (minus the comments)
    elementFileName = pathBaseName + elementFileName;
    this->mElementRawData = this->GetRawDataFromFile (elementFileName);

    // Read the rest of the element data using TokenizeStringsToInts method
    this->mElementData = TokenizeStringsToInts(this->mElementRawData, SPACE_DIM + 1);
    this->mpElementIterator = this->mElementData.begin();

    /*
     * Open edge file and store the lines as a vector of strings (minus the comments)
     * We store edges as "faces" but the superclass
     * provides a GetNextEdgeData method which queries this data.
     */

    edgeFileName = pathBaseName + edgeFileName;
    this->mFaceRawData = this->GetRawDataFromFile (edgeFileName);

    // Read the rest of the face/edge data using TokenizeStringsToInts method
    this->mFaceData = TokenizeStringsToInts(this->mFaceRawData, SPACE_DIM);
    this->mpFaceIterator = this->mFaceData.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::~FemlabMeshReader()
{}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector < std::vector<double> >
    FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToDoubles(std::vector<std::string> rawData)
{
    std::vector < std::vector < double > >tokenized_data; // Output

    // Iterate over the lines of input
    unsigned dimension_count = 0;
    std::vector < std::string >::iterator the_iterator;
    for (the_iterator = rawData.begin(); the_iterator != rawData.end();
         the_iterator++)
    {
        std::string line_of_data = *the_iterator;
        std::stringstream line_stream (line_of_data);

        if (dimension_count == 0)
        {
            // First iteration, build the tokenized_data vector and push in x coordinates
            while (!line_stream.eof())
            {
                double item_coord;

                std::vector < double >x_coord;
                line_stream >> item_coord;
                x_coord.push_back (item_coord);
                tokenized_data.push_back (x_coord);
            }
        }
        else
        {
            unsigned current_node = 0;

            // Other iterations, push in coordinates other than x
            while (!line_stream.eof())
            {
                double item_coord;
                line_stream >> item_coord;
                tokenized_data[current_node].push_back (item_coord);
                current_node++;
            }
        }
        // Dimension of mesh is the same as the line of rawData
        dimension_count++;
    }

    if (SPACE_DIM != dimension_count)
    {
        EXCEPTION("SPACE_DIM  != dimension read from file");
    }
    return (tokenized_data);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector < std::vector < unsigned > >
    FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToInts(std::vector<std::string> rawData, unsigned dimensionOfObject)
{
    std::vector < std::vector < unsigned > >tokenized_data;

    // There are dimensionOfObject lines to be read
    for (unsigned i = 0; i < dimensionOfObject; i++)
    {
        std::string line_of_data = rawData[i];
        std::stringstream line_stream (line_of_data);

        if (i == 0)
        {
            // First iteration, build the tokenized_data vector and push in x coordinates
            while (!line_stream.eof())
            {
                double item_index;

                std::vector < unsigned >first_index;
                line_stream >> item_index;
                first_index.push_back ((unsigned) (item_index - 0.5)); // item indices should be minus 1
                tokenized_data.push_back (first_index);
            }
        }
        else
        {
            unsigned current_node = 0;

            // Other iterations, push in coordinates other than x.
            while (!line_stream.eof())
            {
                double item_index;
                line_stream >> item_index;
                tokenized_data[current_node].
                push_back ((unsigned) (item_index - 0.5));
                current_node++;
            }
        }
    }
    return (tokenized_data);
}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class FemlabMeshReader<1,1>;
template class FemlabMeshReader<1,2>;
template class FemlabMeshReader<1,3>;
template class FemlabMeshReader<2,2>;
template class FemlabMeshReader<2,3>;
template class FemlabMeshReader<3,3>;
