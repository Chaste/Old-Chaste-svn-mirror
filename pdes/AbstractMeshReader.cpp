#include "AbstractMeshReader.hpp"


std::vector<std::string> AbstractMeshReader::GetRawDataFromFile(std::string fileName)
{
	std::vector<std::string> RawDataFromFile;

	std::ifstream dataFile(fileName.c_str());
	
	if (!dataFile.is_open())
	{
		throw Exception("Could not open data file "+fileName+" .");
	}
	
	std::string RawLineFromFile;
	getline(dataFile, RawLineFromFile);
	
	while(dataFile){
		int hashLocation=RawLineFromFile.find('#',0);
		if (hashLocation >= 0)
		{
			RawLineFromFile=RawLineFromFile.substr(0,hashLocation);
		}
		int notBlankLocation=RawLineFromFile.find_first_not_of(" \t",0);
		if (notBlankLocation >= 0)
		{
			RawDataFromFile.push_back(RawLineFromFile);	
		}
		getline(dataFile, RawLineFromFile);
	}


	dataFile.close();
	
	return(RawDataFromFile);
}

int AbstractMeshReader::GetMaxNodeIndex()
{
	std::vector<std::vector<int> >::iterator the_iterator;
	
	int max_node_index = -1;
	
	for(the_iterator = mElementData.begin(); the_iterator < mElementData.end(); the_iterator++)
	{
		std::vector<int> indices = *the_iterator;
		
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

int AbstractMeshReader::GetMinNodeIndex()
{
	std::vector<std::vector<int> >::iterator the_iterator;
	
	int min_node_index = 1000000;
	
	for(the_iterator = mElementData.begin(); the_iterator < mElementData.end(); the_iterator++)
	{
		std::vector<int> indices = *the_iterator;
		
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

std::vector<double> AbstractMeshReader::GetNextNode()
{
	if (mpNodeIterator == mNodeData.end())
	{
		throw Exception("All nodes already got");
	}

	std::vector<double> next_node = *mpNodeIterator;
	
	mpNodeIterator++;
		
	return next_node;
}

std::vector<int> AbstractMeshReader::GetNextElement()
{
	if (mpElementIterator == mElementData.end())
	{
		throw Exception("All elements already got");
	}

	std::vector<int> next_element = *mpElementIterator;
	
	mpElementIterator++;
		
	return next_element;
}

std::vector<int> AbstractMeshReader::GetNextFace()
{
	if (mpFaceIterator == mFaceData.end())
	{
		throw Exception("All faces (or edges) already got");
	}

	std::vector<int> next_face = *mpFaceIterator;
	
	mpFaceIterator++;
		
	return next_face;
}

std::vector<int> AbstractMeshReader::GetNextEdge()
{
	return GetNextFace();
}
