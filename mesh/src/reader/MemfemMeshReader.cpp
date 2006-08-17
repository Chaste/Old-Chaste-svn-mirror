/** Implementation file for the MemfemMeshReader class.*/
#ifndef _MEMFEMMESHREADER_CPP_
#define _MEMFEMMESHREADER_CPP_

#include "MemfemMeshReader.hpp"
/**
 * The constructor takes the base name of a set of Memfem
 * mesh files (ie. the path and name of the files without the suffices)
 * and allows the data to be queried.
 * Typical use:
 *    AbstractMeshReader *spMeshReader=new MemfemMeshReader(
 *		                  "pdes/tests/meshdata/Memfem_slab");
 * Also calls the superclass AbstractMeshReader's constructor
 */

template<int ELEMENT_DIM, int SPACE_DIM>
MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::MemfemMeshReader(std::string pathBaseName)
        : AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>()
{

    //Open node file and store the lines as a vector of strings (minus the comments)
    std::string nodeFileName=pathBaseName+".pts";
    this->mNodeRawData=this->GetRawDataFromFile(nodeFileName);
    
    /* Read single line header which is the number of nodes */
    std::stringstream node_header_stream(this->mNodeRawData[0]);
    unsigned int num_nodes;
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
    
    /* Read single line header which is the number of elements	 */
    std::stringstream element_header_stream(this->mElementRawData[0]);
    unsigned int num_elements;
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

template<int ELEMENT_DIM, int SPACE_DIM>
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
            for (int i = 0; i < SPACE_DIM; i++)
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
 * 	NB: Region markers are currently ignored.
 * Return value is a vector where each item is a vector of ints which represents
 * indices of nodes.
 */

template <int ELEMENT_DIM, int SPACE_DIM>
std::vector<std::vector<int> > MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToInts(
    std::vector<std::string> rawData,
    int dimensionOfObject,
    bool readHeader)
{
    std::vector<std::vector<int> > tokenized_data;
    
    std::vector<std::string>::iterator the_iterator;
    for ( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
    {
        std::string line_of_data=*the_iterator;
        std::stringstream line_stream(line_of_data);
        
        
        if ( readHeader == false || the_iterator!=rawData.begin() )
        {
            std::vector<int> current_indices;
            
            for (int i = 0; i < dimensionOfObject; i++)
            {
                int item_index;
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


template <int ELEMENT_DIM, int SPACE_DIM>
MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::~MemfemMeshReader()
{}
#endif //_MEMFEMMESHREADER_CPP_
