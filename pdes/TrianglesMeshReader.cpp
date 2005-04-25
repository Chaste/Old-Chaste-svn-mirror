/** Implementation file for the TrianglesMeshReader class.*/

#include "TrianglesMeshReader.hpp"
/**
 * The constructor takes the base name of a set of Triangles or
 * Tetgen mesh files (ie. the path and name of the files without the suffices)
 * and allows the data to be queried.
 * Typical use:
 *    AbstractMeshReader *spMeshReader=new TrianglesMeshReader(
 *		                  "pdes/tests/meshdata/disk_522_elements");
 * Also calls the superclass AbstractMeshReader's constructor
 */ 
TrianglesMeshReader::TrianglesMeshReader(std::string pathBaseName)
{
	bool indexed_from_zero = false;
	bool already_checked_indexing = false;
	
	int num_attributes, max_marker;
	
	//Copy path and base name of files to private data
	mPathBaseName=pathBaseName;

	//Open node file and store the lines as a vector of strings (minus the comments) 	
	std::string nodeFileName=pathBaseName+".node";
	mNodeRawData=GetRawDataFromFile(nodeFileName);
	
	/* Read single line header as at:
	 * http://tetgen.berlios.de/fformats.node.html
	 * http://www-2.cs.cmu.edu/~quake/triangle.node.html
	 */
	std::stringstream node_header_stream(mNodeRawData[0]);
	node_header_stream >> mNumNodes >> mDimension >> mNumNodeAttributes >> mMaxNodeBdyMarker;
	
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
	std::string elementFileName=pathBaseName+".ele";
	mElementRawData=GetRawDataFromFile(elementFileName);

 	/* Read single line header as at:
	 * http://tetgen.berlios.de/fformats.ele.html
	 * http://www-2.cs.cmu.edu/~quake/triangle.ele.html
	 */
	 std::stringstream element_header_stream(mElementRawData[0]);
	element_header_stream >> mNumElements >> mNumElementNodes >> mNumElementAttributes;
	
	//Only order-1 triangles or tetrahedra are currently supported
	if (mNumElementNodes != mDimension+1)
	{
		throw Exception("Number of nodes per element is not supported");
	}
	
	// Read the rest of the element data using TokenizeStringsToInts method
	mElementData = TokenizeStringsToInts(mElementRawData,mDimension+1);
 	mpElementIterator = mElementData.begin();
 	
 	
 	//Check that the size of the data matches the information in the header
 	if (mNumElements != mElementData.size())
	{
		throw Exception("Number of elements does not match expected number declared in header");
	}

	/*Open face file and store the lines as a vector of strings (minus the comments) 	
	 * Note that the Triangles equivalent of a face file is an edge file.
	 * The two file formats are similar.  We store edges as "faces" but the superclass
	 * provides a GetNextEdge method which queries this data.
	 */

	std::string faceFileName;

	if (mDimension == 2)
	{
		faceFileName=pathBaseName+".edge";
	}
	else if (mDimension == 3)
	{
		faceFileName=pathBaseName+".face";
	}
	
	mFaceRawData=GetRawDataFromFile(faceFileName);
	
	/* Read single line header as at:
	 * http://tetgen.berlios.de/fformats.face.html
	 * http://www-2.cs.cmu.edu/~quake/triangle.edge.html
	 */
	std::stringstream face_header_stream(mFaceRawData[0]);
	face_header_stream >> mNumFaces >> mMaxFaceBdyMarker;

	//mNumBoundaryFaces = mNumFaces; //temporary

	// Read the rest of the face/edge data using TokenizeStringsToInts method
	mFaceData = TokenizeStringsToInts(mFaceRawData,mDimension);
	mpFaceIterator = mFaceData.begin();
	
//	std::stringstream boundary_face_header_stream(mFaceRawData[0]);
//	boundary_face_header_stream >> mNumBoundaryFaces;
	mBoundaryFaceData = CullInternalFaces();
	mpBoundaryFaceIterator = mBoundaryFaceData.begin();
	
	//Check that the size of the data matches the information in the header
	if (mNumFaces != mFaceData.size())
	{
		throw Exception("Number of faces does not match expected number declared in header");
	}
	
	
	
}

/** 
 * TokenizeStringsToDoubles is specific to reading node data which came from 
 * a Triangles or Tetgen file.
 * http://tetgen.berlios.de/fformats.node.html
 * http://www-2.cs.cmu.edu/~quake/triangle.node.html
 *  Each string is expected to be an index number, 2 or 3 doubles (representing x,y,z)
 *  an optional list of attributes and an optional boundary marker
 * 	NB: Attributes and boundary markers are currently ignored.
 * Return value is a vector where each item is a vector of double which represents 
 * position.  Indices are implicit in the vector.
 */
	
std::vector<std::vector<double> > TrianglesMeshReader::TokenizeStringsToDoubles(
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
     		int item_number;
     		
     		//Read item index     		
     		line_stream >> item_number;
     		
     		//Watch for item zero (this will only happen at most once!)
     		if (item_number == 0)
     		{
     			mIndexFromZero = true;
     		}
     		
     		//Throw an error if the item numbers are out of order
     		if (mIndexFromZero == true)
     		{
     			if (item_number != tokenized_data.size())
     			{
     				std::stringstream item_number_as_string;
					item_number_as_string << item_number;
     				throw Exception("Node number " + item_number_as_string.str() + " is out of order in input file.");
     			}
     		}
     		else
     		{
     			if (item_number-1 != tokenized_data.size())
     			{
     				std::stringstream item_number_as_string;
					item_number_as_string << item_number;
     				throw Exception("Node number " + item_number_as_string.str() + " is out of order in input file.");
     			}
     		}
     		
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
 * 	NB: Attributes and boundary markers are currently ignored.
 * Return value is a vector where each item is a vector of ints which represents 
 * indices of nodes.  
 */
std::vector<std::vector<int> > TrianglesMeshReader::TokenizeStringsToInts(
												std::vector<std::string> rawData,
												int dimensionOfObject)
{
	std::vector<std::vector<int> > tokenized_data;
	
	std::vector<std::string>::iterator the_iterator;
   	for( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
   	{
     	std::string line_of_data=*the_iterator;
     	std::stringstream line_stream(line_of_data);
     	
     	if (the_iterator!=rawData.begin())
     	{
     		std::vector<int> current_indices;
     		int item_number;
     		
     		line_stream >> item_number;
     		
     		for (int i = 0; i < dimensionOfObject; i++)
     		{
     			int item_index; 
     			line_stream >> item_index;
     			//If the nodes have been indexed from one then we need to shift the indices
     			if (mIndexFromZero == false)
     			{
     				item_index -= 1;
     			}
     			current_indices.push_back(item_index);
     		}
     		
     		tokenized_data.push_back(current_indices);
     	}
     	
    }
    
    return tokenized_data;
}


/**
 * Destructor
 */
TrianglesMeshReader::~TrianglesMeshReader()
{
}
