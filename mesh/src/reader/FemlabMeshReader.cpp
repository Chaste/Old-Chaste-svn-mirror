/** Implementation file for the FemlabMeshReader class.*/
#ifndef _FEMLABMESHREADER_CPP_
#define _FEMLABMESHREADER_CPP_

#include "FemlabMeshReader.hpp"
#include "Exception.hpp"

/**
 * The constructor takes the path to and names of a set of Femlab mesh files
 * (ie. the node, elements and face files (in that order) and allows the data to
 * be queried.
 * Typical use:
 *    AbstractMeshReader *spMeshReader=new FemlabMeshReader(
 * 						  "pdes/tests/meshdata/",
 *		                  "femlab_lshape_nodes.dat",
 * 						  "femlab_lshape_elements.dat",
 * 						  "femlab_lshape_edges.dat",);
 * Also calls the superclass AbstractMeshReader's constructor
 */
template <int ELEMENT_DIM, int SPACE_DIM>
FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::FemlabMeshReader (std::string pathBaseName,
                                                            std::string nodeFileName,
                                                            std::string elementFileName,
                                                            std::string edgeFileName)
{

    //Open node file and store the lines as a vector of strings (minus the comments)
    nodeFileName = pathBaseName + nodeFileName;
    this->mNodeRawData = this->GetRawDataFromFile (nodeFileName);
    
    // Read the node data using TokenizeStringsToDoubles method
    this->mNodeData = TokenizeStringsToDoubles (this->mNodeRawData);
    
    //Initialise iterator for public GetNextNode method
    this->mpNodeIterator = this->mNodeData.begin ();
    
    
    //Open element file and store the lines as a vector of strings (minus the comments)
    elementFileName = pathBaseName + elementFileName;
    this->mElementRawData = this->GetRawDataFromFile (elementFileName);
    
    // Read the rest of the element data using TokenizeStringsToInts method
    this->mElementData = TokenizeStringsToInts (this->mElementRawData, SPACE_DIM + 1);
    this->mpElementIterator = this->mElementData.begin ();
    
    
    /*Open edge file and store the lines as a vector of strings (minus the comments)
     * We store edges as "faces" but the superclass
     * provides a GetNextEdge method which queries this data.
     */
    
    edgeFileName = pathBaseName + edgeFileName;
    this->mFaceRawData = this->GetRawDataFromFile (edgeFileName);
    
    // Read the rest of the face/edge data using TokenizeStringsToInts method
    this->mFaceData = TokenizeStringsToInts (this->mFaceRawData, SPACE_DIM);
    this->mpFaceIterator = this->mFaceData.begin ();
}

/**
 * TokenizeStringsToDoubles is specific to reading node data which came from
 * a Femlab or Matlab PDE toolbox file.
 *
 * Each string is expected to be a series of doubles.
 * Return value is a vector where each item is a vector of double which represents
 * position.  Indices are implicit in the vector.
 */

template <int ELEMENT_DIM, int SPACE_DIM>
std::vector < std::vector < double > >
FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToDoubles (std::vector < std::string >
        rawData)
{
    std::vector < std::vector < double > >tokenized_data;        // Output
    
    //Iterate over the lines of input
    unsigned dimension_count = 0;
    std::vector < std::string >::iterator the_iterator;
    for (the_iterator = rawData.begin (); the_iterator != rawData.end ();
         the_iterator++)
    {
        std::string line_of_data = *the_iterator;
        std::stringstream line_stream (line_of_data);
        
        if (dimension_count == 0)
        {
            //First iteration, build the tokenized_data vector and push in x coordinates
            while (!line_stream.eof ())
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
            
            //Other iterations, push in coordinates other than x.
            while (!line_stream.eof ())
            {
                double item_coord;
                line_stream >> item_coord;
                tokenized_data[current_node].push_back (item_coord);
                current_node++;
            }
        }
        //dimension of mesh is the same as the line of rawData.
        dimension_count++;
        
    }
    
    if (SPACE_DIM != dimension_count)
    {
        EXCEPTION("SPACE_DIM  != dimension read from file");
    }
    return (tokenized_data);
}


/**
 * TokenizeStringsToInts is for reading element, face or edge data which came from
 * a Femlab or Matlab PDE toolbox file.
 *  Each string is expected to be a series of unsigned which represent:
 *  The first several lines denote the indices of nodes
 *  The rest contains extra information which are ignored currently.
 *  ( In 2-D: 2 indices for an edge, 3 for a triangle)
 *  ( In 3-D: 3 indices for a face, 4 for a tetrahedron)
 * Return value is a vector where each item is a vector of ints which represents
 * indices of nodes.
 */
template <int ELEMENT_DIM, int SPACE_DIM>
std::vector < std::vector < unsigned > >
FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToInts (std::vector < std::string > rawData,
        unsigned dimensionOfObject)
{
    std::vector < std::vector < unsigned > >tokenized_data;
    
    /* There are dimensionOfObject lines to be read */
    for (unsigned i = 0; i < dimensionOfObject; i++)
    {
        std::string line_of_data = rawData[i];
        std::stringstream line_stream (line_of_data);
        
        if (i == 0)
        {
            //First iteration, build the tokenized_data vector and push in x coordinates
            while (!line_stream.eof ())
            {
                double item_index;
                
                std::vector < unsigned >first_index;
                line_stream >> item_index;
                first_index.push_back ((unsigned) (item_index - 0.5));       //item indices should be minus 1.
                tokenized_data.push_back (first_index);
            }
        }
        else
        {
            unsigned current_node = 0;
            
            //Other iterations, push in coordinates other than x.
            while (!line_stream.eof ())
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


/**
 * Destructor
 */
template <int ELEMENT_DIM, int SPACE_DIM>
FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::~FemlabMeshReader ()
{}
#endif //_FEMLABMESHREADER_CPP_
