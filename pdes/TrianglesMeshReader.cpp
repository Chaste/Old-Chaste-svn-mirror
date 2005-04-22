#include "TrianglesMeshReader.hpp"

TrianglesMeshReader::TrianglesMeshReader(std::string pathBaseName)
{
	bool indexed_from_zero = false;
	bool already_checked_indexing = false;
	
	int num_attributes, max_marker;
	
	mPathBaseName=pathBaseName;

	std::string nodeFileName=pathBaseName+".node";
	mNodeRawData=GetRawDataFromFile(nodeFileName);
	
	std::stringstream node_header_stream(mNodeRawData[0]);
	node_header_stream >> mNumNodes >> mDimension >> mNumNodeAttributes >> mMaxNodeBdyMarker;
	
	mNodeData = TokenizeStringsToDoubles(mNodeRawData);
	mpNodeIterator = mNodeData.begin();
	
	if (mNumNodes != mNodeData.size())
	{
		throw Exception("Number of nodes does not match expected number declared in header");
	}
	
	
	std::string elementFileName=pathBaseName+".ele";
	mElementRawData=GetRawDataFromFile(elementFileName);

 	std::stringstream element_header_stream(mElementRawData[0]);
	element_header_stream >> mNumElements >> mNumElementNodes >> mNumElementAttributes;
	
	if (mNumElementNodes != mDimension+1)
	{
		throw Exception("Number of nodes per element is not supported");
	}
	
	mElementData = TokenizeStringsToInts(mElementRawData,mDimension+1);
 	mpElementIterator = mElementData.begin();
 	
 	
 	
 	if (mNumElements != mElementData.size())
	{
		throw Exception("Number of elements does not match expected number declared in header");
	}

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
	
	std::stringstream face_header_stream(mFaceRawData[0]);
	face_header_stream >> mNumFaces >> mMaxFaceBdyMarker;
	
	mFaceData = TokenizeStringsToInts(mFaceRawData,mDimension);
	mpFaceIterator = mFaceData.begin();
	
	if (mNumFaces != mFaceData.size())
	{
		throw Exception("Number of faces does not match expected number declared in header");
	}
	
	
	
}
	
std::vector<std::vector<double> > TrianglesMeshReader::TokenizeStringsToDoubles(
												std::vector<std::string> rawData)
{
	std::vector<std::vector<double> > tokenized_data;
	
	std::vector<std::string>::iterator the_iterator;
   	for( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
   	{
     	std::string line_of_data=*the_iterator;
     	std::stringstream line_stream(line_of_data);
     	
     	if (the_iterator!=rawData.begin())
     	{
     		std::vector<double> current_coords;
     		int item_number;
     		     		
     		line_stream >> item_number;
     		
     		if (item_number == 0)
     		{
     			mIndexFromZero = true;
     		}
     		
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
     		
     		for (int i = 0; i < mDimension; i++)
     		{
     			double item_coord; 
     			line_stream >> item_coord;
     			current_coords.push_back(item_coord);
     		}
     		
     		tokenized_data.push_back(current_coords);
     	}
     	
    }
    
    return tokenized_data;
}


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



TrianglesMeshReader::~TrianglesMeshReader()
{
}
