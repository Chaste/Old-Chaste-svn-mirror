#include "TrianglesMeshReader.hpp"

TrianglesMeshReader::TrianglesMeshReader(std::string pathBaseName)
{
	
	mPathBaseName=pathBaseName;

	std::string nodeFileName=pathBaseName+".node";
	mNodeRawData=GetRawDataFromFile(nodeFileName);

	std::vector<std::string>::iterator theNodeIterator;
   	for( theNodeIterator = mNodeRawData.begin(); theNodeIterator != mNodeRawData.end(); theNodeIterator++ )
   	{
     	std::string lineOfData=*theNodeIterator;
     	std::stringstream lineStream(lineOfData);
     	
     	if (theNodeIterator==mNodeRawData.begin())
     	{
     		lineStream >> mNumNodes;     
     	}
     	else
     	{
     		std::vector<double> currentNodeCoords;
     		int nodeNumber, nodeCoords[3]; 
     		
     		lineStream >> nodeNumber >> nodeCoords[0] >> nodeCoords[1] >> nodeCoords[2];
     		
     		for (int i = 0; i < 3; i++)
     		{
     			currentNodeCoords.push_back(nodeCoords[i]);
     		}
     		
     		mNodeData.push_back(currentNodeCoords);
     	}
     	
    }
    
    if (mNumNodes != mNodeData.size())
	{
		throw Exception("Number of nodes does not match expected number declared in header");
	}





	std::string elementFileName=pathBaseName+".ele";
	mElementRawData=GetRawDataFromFile(elementFileName);

 	std::vector<std::string>::iterator theElementIterator;
   	for( theElementIterator = mElementRawData.begin(); theElementIterator != mElementRawData.end(); theElementIterator++ )
   	{
     	std::string lineOfData=*theElementIterator;
     	std::stringstream lineStream(lineOfData);
     	
     	if (theElementIterator==mElementRawData.begin())
     	{
     		lineStream >> mNumElements;     
     	}
     	else
     	{
     		std::vector<int> currentElementNodes;
     		int elementNumber, elementNodes[3]; 
     		
     		lineStream >> elementNumber >> elementNodes[0] >> elementNodes[1] >> elementNodes[2];
     		
     		for (int i = 0; i < 3; i++)
     		{
     			currentElementNodes.push_back(elementNodes[i]);
     		}
     		
     		mElementData.push_back(currentElementNodes);
     	}
     	
    }
    
    if (mNumElements != mElementData.size())
	{
		throw Exception("Number of elements does not match expected number declared in header");
	}



	std::string faceFileName=pathBaseName+".edge";
	mFaceRawData=GetRawDataFromFile(faceFileName);
	
	std::vector<std::string>::iterator theFaceIterator;
   	for( theFaceIterator = mFaceRawData.begin(); theFaceIterator != mFaceRawData.end(); theFaceIterator++ )
   	{
     	std::string lineOfData=*theFaceIterator;
     	std::stringstream lineStream(lineOfData);
     	
     	if (theFaceIterator==mFaceRawData.begin())
     	{
     		lineStream >> mNumFaces;     
     	}
     	else
     	{
     		std::vector<int> currentFaceNodes;
     		int faceNumber, faceNodes[2]; 
     		
     		lineStream >> faceNumber >> faceNodes[0] >> faceNodes[1];
     		
     		for (int i = 0; i < 2; i++)
     		{
     			currentFaceNodes.push_back(faceNodes[i]);
     		}
     		
     		mFaceData.push_back(currentFaceNodes);
     	}
     	
    }
        
    if (mNumFaces != mFaceData.size())
	{
		throw Exception("Number of faces does not match expected number declared in header");
	}





		
	//int numAttributes, maxBoundaryMarkers;
	//nodeFile >> mNumNodes >> mDimension >> numAttributes >> maxBoundaryMarkers;
	
}

TrianglesMeshReader::~TrianglesMeshReader()
{
}
