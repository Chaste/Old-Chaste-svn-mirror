/** Implementation file for the MemfemMeshReader class.*/

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


MemfemMeshReader::MemfemMeshReader(std::string pathBaseName)
{
	
	//Copy path and base name of files to private data
	mPathBaseName=pathBaseName;

	//Open node file and store the lines as a vector of strings (minus the comments) 	
	std::string nodeFileName=pathBaseName+".pts";
	mNodeRawData=GetRawDataFromFile(nodeFileName);
	
	/* Read single line header which is the number of nodes */
	std::stringstream node_header_stream(mNodeRawData[0]);
	node_header_stream >> mNumNodes;
	
	/* All Memfem data is in 3-d. */
	mDimension = 3; 
	
	// Read the rest of the node data using TokenizeStringsToDoubles method
	mNodeData = TokenizeStringsToDoubles(mNodeRawData);
	//Initialise iterator for public GetNextNode method
	mpNodeIterator = mNodeData.begin();
	
	//Check that the size of the data matches the information in the header
	
	if (mNumNodes != mNodeData.size())
	{
		throw Exception("Number of nodes does not match expected number declared in header");
	}
	
	//Open element file and store the lines as a vector of strings (minus the comments) 	
	std::string elementFileName=pathBaseName+".tetras";
	mElementRawData=GetRawDataFromFile(elementFileName);

 	/* Read single line header which is the number of elements	 */
	 std::stringstream element_header_stream(mElementRawData[0]);
	element_header_stream >> mNumElements;
	

	// Read the rest of the element data using TokenizeStringsToInts method
	mElementData = TokenizeStringsToInts(mElementRawData,mDimension+1, true);
 	mpElementIterator = mElementData.begin();
 	
 	
 	//Check that the size of the data matches the information in the header
  	if (mNumElements != mElementData.size())
	{
		throw Exception("Number of elements does not match expected number declared in header");
	}

	//Open boundary face file and store the lines as a vector of strings (minus the comments) 	
	std::string faceFileName=pathBaseName+".tri";
	mFaceRawData=GetRawDataFromFile(faceFileName);
	
	/* There is no header.
	 */
	
	// Read the face/edge data using TokenizeStringsToInts method
	mBoundaryFaceData = TokenizeStringsToInts(mFaceRawData,mDimension,false);
	mpBoundaryFaceIterator = mBoundaryFaceData.begin();
	mNumBoundaryFaces = mBoundaryFaceData.size();
	
	
	
}


/** 
 * TokenizeStringsToDoubles is specific to reading node data which came from 
 * a Memfem file.
 * Each string is expected to be 3 doubles (representing x,y,z)
 * Return value is a vector where each item is a vector of doubles which represents 
 * position.  Indices are implicit in the vector.
 */
	
std::vector<std::vector<double> > MemfemMeshReader::TokenizeStringsToDoubles(
												std::vector<std::string> rawData)
{
	std::vector<std::vector<double> > tokenized_data; // Output
	
	//Iterate over the lines of input
	std::vector<std::string>::iterator the_iterator;
   	for( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
   	{
     	std::string line_of_data=*the_iterator;
     	std::stringstream line_stream(line_of_data);
     	
     	if (the_iterator!=rawData.begin()) //Ignore the header string
     	{
     		std::vector<double> current_coords;
     		
     		//Form the vector which represents the position of this item
     		for (int i = 0; i < mDimension; i++)
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
std::vector<std::vector<int> > MemfemMeshReader::TokenizeStringsToInts(
												std::vector<std::string> rawData,
												int dimensionOfObject,
												bool readHeader)
{
	std::vector<std::vector<int> > tokenized_data;
	
	std::vector<std::string>::iterator the_iterator;
   	for( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
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



MemfemMeshReader::~MemfemMeshReader()
{
}
