#include "AbstractMeshReader.hpp"
#include "common/Exception.hpp"

// AbstractMeshReader.cpp

 
/**
 * Reads an input file fileName, removes comments (indicated by a #) and blank 
 * lines and returns a vector of strings. Each string corresponds to one line
 * of the input file.
 * 
 * 
 */

std::vector<std::string> AbstractMeshReader::GetRawDataFromFile(std::string fileName)
{
	// Open raw data file
	
	std::vector<std::string> RawDataFromFile;
	std::ifstream dataFile(fileName.c_str());
	
	
	
	/**
	 * Checks that input file has been opened correctly. If not throws an 
	 * exception that should be caught by the user.
	 * 
	 */
	
	if (!dataFile.is_open())
	{
		throw Exception("Could not open data file "+fileName+" .");
	}
	
	
	
	// Read each line in turn
	
	std::string RawLineFromFile;
	getline(dataFile, RawLineFromFile);
	
	while(dataFile){
		
		// Remove comments
		
		int hashLocation=RawLineFromFile.find('#',0);
		if (hashLocation >= 0)
		{
			RawLineFromFile=RawLineFromFile.substr(0,hashLocation);
		}
		
	
		// Remove blank lines
		
		int notBlankLocation=RawLineFromFile.find_first_not_of(" \t",0);
		if (notBlankLocation >= 0)
		{
			RawDataFromFile.push_back(RawLineFromFile);	
		}
		
		
		// Move onto next line
		
		getline(dataFile, RawLineFromFile);
	}


	dataFile.close(); // Closes the data file
	
	return(RawDataFromFile);
}



/**
 * Returns the maximum node index. Used in testing to check that output nodes 
 * are always indexed from zero even if they are input indexed from one.
 * 
 */

int AbstractMeshReader::GetMaxNodeIndex()
{
	//Initialize an interator for the vector of nodes
	std::vector<std::vector<int> >::iterator the_iterator;
	
	int max_node_index = -1; // Must be negative
	
	for(the_iterator = mElementData.begin(); the_iterator < mElementData.end(); the_iterator++)
	{
		std::vector<int> indices = *the_iterator; // the_iterator points at each line in turn
		
		for(int i = 0; i < mDimension+1; i++)
		{
			if(indices[i] > max_node_index)
			{
				max_node_index = indices[i];
			}
		}
	}
	
	return max_node_index;
}



/**
 * Returns the minimum node index. Used in testing to check that output nodes 
 * are always indexed from zero even if they are input indexed from one.
 * 
 */

int AbstractMeshReader::GetMinNodeIndex()
{
	//Initialize an interator for the vector of nodes
	std::vector<std::vector<int> >::iterator the_iterator;
	
	int min_node_index = 1000000; // A large integer
	
	for(the_iterator = mElementData.begin(); the_iterator < mElementData.end(); the_iterator++)
	{
		std::vector<int> indices = *the_iterator; // the_iterator points at each line in turn
		
		for(int i = 0; i < mDimension+1; i++)
		{
			if(indices[i] < min_node_index)
			{
				min_node_index = indices[i];
			}
		}
	}
	
	return min_node_index;
}



/**
 * Returns a vector of the coordinates of each node in turn, starting with
 * node 0 the first time it is called followed by nodes 1, 2, ... , mNumNodes-1.
 * 
 */

std::vector<double> AbstractMeshReader::GetNextNode()
{
	/**
	 * Checks that there are still some nodes left to read. If not throws an
	 * exception that must be caught by the user.
	 * 
	 */
	
	if (mpNodeIterator == mNodeData.end())
	{
		throw Exception("All nodes already got");
	}

	std::vector<double> next_node = *mpNodeIterator;
	
	mpNodeIterator++;
		
	return next_node;
}



/**
 * Returns a vector of the nodes of each element in turn, starting with
 * element 0 the first time it is called followed by elements 1, 2, ... , 
 * mNumElements-1.
 * 
 */

std::vector<int> AbstractMeshReader::GetNextElement()
{
	/**
	 * Checks that there are still some elements left to read. If not throws an
	 * exception that must be caught by the user.
	 * 
	 */
	
	if (mpElementIterator == mElementData.end())
	{
		throw Exception("All elements already got");
	}

	std::vector<int> next_element = *mpElementIterator;
	
	mpElementIterator++;
		
	return next_element;
}

/**
 * 
 */

void AbstractMeshReader::Reset()
{	
	mpElementIterator = mElementData.begin();
    mpFaceIterator = mFaceData.begin();
    mpBoundaryFaceIterator = mBoundaryFaceData.begin();
    mpNodeIterator = mNodeData.begin();
}


/**
 * Returns a vector of the nodes of each face in turn, starting with face 0 the 
 * first time it is called followed by faces 1, 2, ... , mNumFaces-1.
 * 
 * Is a synonum of GetNextEdge(). The two functions can be used interchangeably,
 * i.e. they use the same iterator.
 * 
 */

std::vector<int> AbstractMeshReader::GetNextFace()
{
	/**
	 * Checks that there are still some faces left to read. If not throws an
	 * exception that must be caught by the user.
	 * 
	 */
	
	if (mpFaceIterator == mFaceData.end())
	{
		throw Exception("All faces (or edges) already got");
	}

	std::vector<int> next_face = *mpFaceIterator;
	
	mpFaceIterator++;
		
	return next_face;
}

std::vector<int> AbstractMeshReader::GetNextBoundaryFace()
{
	/**
	 * Checks that there are still some faces left to read. If not throws an
	 * exception that must be caught by the user.
	 * 
	 */
	
	if (mpBoundaryFaceIterator == mBoundaryFaceData.end())
	{
		throw Exception("All faces (or edges) already got");
	}

	std::vector<int> next_face = *mpBoundaryFaceIterator;
	
	mpBoundaryFaceIterator++;
		
	return next_face;
}



/**
 * Returns a vector of the nodes of each edge in turn, starting with edge 0 the 
 * first time it is called followed by edges 1, 2, ... , mNumFaces-1.
 * 
 * Is a synonum of GetNextFace(). The two functions can be used interchangeably,
 * i.e. they use the same iterator.
 * 
 */

std::vector<int> AbstractMeshReader::GetNextEdge()
{
	// Call GetNextFace()
	return GetNextFace();
}

std::vector<int> AbstractMeshReader::GetNextBoundaryEdge()
{
	// Call GetNextFace()
	return GetNextBoundaryFace();
}

/**
 * Some file formats have a list of all faces and others only have a
 * list of the boundary faces. In order for consistency among the 
 * concrete classes, we here provide the functionality to remove the 
 * internal faces.
 * The original face data is actually preserved. We simply copy the 
 * boundary face data to another place.
 */

std::vector< std::vector<int> > AbstractMeshReader::CullInternalFaces()
{
	std::vector< std::vector<int> > boundary_faces; 	

	
	// Iterate over all faces
	int num_faces=GetNumFaces();
	for (int faceIndex=0; faceIndex<num_faces; faceIndex++)
	{
		std::vector<int> current_face = mFaceData[faceIndex];
		// Initialise count of how many elements this face belongs to
		int num_of_owning_elements = 0;

		// Iterate over all elements
		int num_elements=GetNumElements();		
		for (int elementIndex = 0; elementIndex < num_elements && num_of_owning_elements < 2; 
															elementIndex++ )
		{
			int num_of_matches = 0;
			std::vector<int> current_element = mElementData[elementIndex];
			
			//Iterate over indices of the element
			for (int count = 0; count<=mDimension && num_of_matches >= count-1; count++)
			{
				//Iterate over indices of the face
				for (int i=0; i<mDimension; i++)
				{
					if (current_face[i] == current_element[count])
					{
						num_of_matches++;
						break;
					}
				}						
			}	
			//The current face is a member of the current element
			if (num_of_matches == mDimension)
			{
				num_of_owning_elements++;
			}		

		}
		// A face belonging to exactly one element is a boundary face
		if ( num_of_owning_elements == 1 )
		{
			boundary_faces.push_back(current_face);
			//mNumBoundaryFaces++;
		}		
		else if (num_of_owning_elements != 2 )
		{
			throw Exception("All faces should belong to either one or two elements. ");
		}
	}	
	
	return boundary_faces;
}
