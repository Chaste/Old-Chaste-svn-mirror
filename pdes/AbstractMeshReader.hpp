#ifndef _ABSTRACTMESHREADER_HPP_
#define _ABSTRACTMESHREADER_HPP_

#include <string>
#include <fstream>
#include <iostream>

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
		
		double **mppNodes;
		int **mppElements;
	
		int** MakeIntDynamicArray(int numRows, int numColumns)
		{
			int **ppDynamicArray = new int *[numRows];
			for (int i = 0; i < numRows; i++)
			{
				ppDynamicArray[i] = new int [numColumns];
				for (int j = 0; j < numColumns; j++)
				{
					ppDynamicArray[i][j] = 0;	
				}
			}
			return ppDynamicArray;
		}
		
		double** MakeDoubleDynamicArray(int numRows, int numColumns)
		{
			double **ppDynamicArray = new double *[numRows];
			for (int i = 0; i < numRows; i++)
			{
				ppDynamicArray[i] = new double [numColumns];
				for (int j = 0; j < numColumns; j++)
				{
					ppDynamicArray[i][j] = 0.0;	
				}
			}
			return ppDynamicArray;
		}
	
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
