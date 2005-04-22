#ifndef _ABSTRACTMESHREADER_HPP_
#define _ABSTRACTMESHREADER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Exception.hpp"

class AbstractMeshReader
{
	//private:
	protected:
		int mNumElements;
		int mNumNodes;
		int mNumFaces;
		int mDimension;
		
		int mNumNodeAttributes;
		int mMaxNodeBdyMarker;
		int mNumElementNodes;
		int mNumElementAttributes;
		int mMaxFaceBdyMarker;
		
		std::string mPathBaseName;
		std::vector<std::string> mNodeRawData;  //Contents of file with comments removed		
		std::vector<std::string> mElementRawData;  //Contents of file with comments removed		
		std::vector<std::string> mFaceRawData;  //Contents of file with comments removed		
	
		std::vector< std::vector<double> > mNodeData;
		std::vector< std::vector<int> > mElementData;
		std::vector< std::vector<int> > mFaceData;
		
		std::vector< std::vector<double> >::iterator mpNodeIterator;
		std::vector< std::vector<int> >::iterator mpElementIterator;
		std::vector< std::vector<int> >::iterator mpFaceIterator;
	
		bool mIndexFromZero;
	
		std::vector<std::string> GetRawDataFromFile(std::string fileName);
		
	public:
		AbstractMeshReader()
		{
			mNumElements = 0;
			mNumNodes = 0;
			mNumFaces = 0;
			mDimension = 0;
			
			mNumNodeAttributes = 0;
			mMaxNodeBdyMarker = 0;
			mNumElementNodes = 0;
			mNumElementAttributes = 0;
			mMaxFaceBdyMarker = 0;
			
			mPathBaseName = "";
			mIndexFromZero = false;
		}
		int GetNumElements(){return mNumElements;}
		int GetNumNodes(){return mNumNodes;}
		int GetNumFaces(){return mNumFaces;}
		int GetNumEdges(){return mNumFaces;}		
		int GetDimension(){return mDimension;}		
		
		int GetMaxNodeIndex();
		int GetMinNodeIndex();
		
		std::vector<double> GetNextNode();
		std::vector<int> GetNextElement();
		std::vector<int> GetNextEdge();
		std::vector<int> GetNextFace();

};
#endif //_ABSTRACTMESHREADER_HPP_
