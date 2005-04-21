#include "TrianglesMeshReader.hpp"

TrianglesMeshReader::TrianglesMeshReader(std::string pathBaseName)
{
	
	mPathBaseName=pathBaseName;

	std::string nodeFileName=pathBaseName+".node";
	mNodeRawData=GetRawDataFromFile(nodeFileName);







	std::string elementFileName=pathBaseName+".ele";
	mElementRawData=GetRawDataFromFile(elementFileName);

 	std::vector<std::string>::iterator theIterator;
   	for( theIterator = mElementRawData.begin(); theIterator != mElementRawData.end(); theIterator++ )
   	{
     	std::string lineOfData=*theIterator;
     	std::stringstream lineStream(lineOfData);
     	
     	if (theIterator==mElementRawData.begin())
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

	
	//int numAttributes, maxBoundaryMarkers;
	//nodeFile >> mNumNodes >> mDimension >> numAttributes >> maxBoundaryMarkers;
	
}

TrianglesMeshReader::~TrianglesMeshReader()
{
}
