
#ifndef _TRIANGLESMESHREADER_CPP_
#define _TRIANGLESMESHREADER_CPP_

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
template <int ELEMENT_DIM, int SPACE_DIM>
TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::TrianglesMeshReader(std::string pathBaseName)
{	
	//Open node file and store the lines as a vector of strings (minus the comments) 	
	std::string nodeFileName=pathBaseName+".node";
	this->mNodeRawData=this->GetRawDataFromFile(nodeFileName);
	
	/* Read single line header as at:
	 * http://tetgen.berlios.de/fformats.node.html
	 * http://www-2.cs.cmu.edu/~quake/triangle.node.html
	 */
	std::stringstream node_header_stream(this->mNodeRawData[0]);
	unsigned int num_nodes;
    unsigned int dimension;
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
		EXCEPTION("Number of nodes does not match expected number declared in header");
	}
 
    if (ELEMENT_DIM < SPACE_DIM)
    {
        
        ReadFacesAsElements(pathBaseName);
        
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
	unsigned int num_elements;
	element_header_stream >> num_elements >> this->mNumElementNodes >> this->mNumElementAttributes;
	
	//Only order 1 triangles or tetrahedra are currently supported
	if (this->mNumElementNodes != ELEMENT_DIM+1)
	{
		EXCEPTION("Number of nodes per element is not supported");
	}

	// Read the rest of the element data using TokenizeStringsToInts method
	this->mElementData = TokenizeStringsToInts(this->mElementRawData,ELEMENT_DIM+1);
 	this->mpElementIterator = this->mElementData.begin();
 	
 	//Check that the size of the data matches the information in the header
 	if (num_elements != this->mElementData.size())
	{
		EXCEPTION("Number of elements does not match expected number declared in header");
	}

	/*Open face file and store the lines as a vector of strings (minus the comments) 	
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
        int num_faces = this->GetNumNodes();
        for (int i=0; i<num_faces; i++)
        {
            std::vector<int> current_item;
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
	unsigned int num_faces;
	face_header_stream >> num_faces >> this->mMaxFaceBdyMarker;

	//mNumBoundaryFaces = mNumFaces; //temporary

	// Read the rest of the face/edge data using TokenizeStringsToInts method
	this->mFaceData = TokenizeStringsToInts(this->mFaceRawData,ELEMENT_DIM);
	this->mpFaceIterator = this->mFaceData.begin();
	
	//Check that the size of the data matches the information in the header
	if (num_faces != this->mFaceData.size())
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
 * 	NB: Attributes and boundary markers are currently ignored.
 * Return value is a vector where each item is a vector of doubles which represents 
 * position.  Indices are implicit in the vector.
 */
	
template <int ELEMENT_DIM, int SPACE_DIM>
std::vector<std::vector<double> > TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToDoubles(
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
     		unsigned int item_number;
     		
     		//Read item index     		
     		line_stream >> item_number;
     		
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
     				std::stringstream item_number_as_string;
					item_number_as_string << item_number;
     				EXCEPTION("Node number " + item_number_as_string.str() + " is out of order in input file.");
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
template <int ELEMENT_DIM, int SPACE_DIM>
std::vector<std::vector<int> > TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToInts(
												std::vector<std::string> rawData,
												int dimensionOfObject)
{
	std::vector< std::vector<int> > tokenized_data;
	
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
     			if (this->mIndexFromZero == false)
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

template <int ELEMENT_DIM, int SPACE_DIM>
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
    unsigned int num_elements;
    face_header_stream >> num_elements >> this->mMaxFaceBdyMarker;
    
    // Read the rest of the element data using TokenizeStringsToInts method
    this->mElementData = TokenizeStringsToInts(this->mElementRawData,ELEMENT_DIM+1);
    this->mpElementIterator = this->mElementData.begin();
    //Check that the size of the data matches the information in the header
    
    if (num_elements != this->mElementData.size())
    {
        EXCEPTION("Number of elements does not match expected number declared in header");
    }
}

/**
 * Destructor
 */
template <int ELEMENT_DIM, int SPACE_DIM>
TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::~TrianglesMeshReader()
{
}

#endif // _TRIANGLESMESHREADER_CPP_
