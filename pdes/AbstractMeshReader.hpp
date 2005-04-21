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
		std::string mPathBaseName;
		std::vector<std::string> mNodeRawData;  //Contents of file with comments removed		
		std::vector<std::string> mElementRawData;  //Contents of file with comments removed		
		std::vector<std::string> mFaceRawData;  //Contents of file with comments removed		
	
		std::vector< std::vector<int> > mElementData;
	
		std::vector<std::string> GetRawDataFromFile(std::string fileName);
		
	public:
		AbstractMeshReader()
		{
			mNumElements = 0;
			mNumNodes = 0;
			mNumFaces = 0;
			mDimension = 0;
			mPathBaseName = "";
		}
		int GetNumElements(){return mNumElements;}
		int GetNumNodes(){return mNumNodes;}
		int GetNumFaces(){return mNumFaces;}
		int GetDimension(){return mDimension;}		

};
#endif //_ABSTRACTMESHREADER_HPP_
