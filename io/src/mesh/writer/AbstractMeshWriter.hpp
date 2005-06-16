#ifndef _ABSTRACTMESHWRITER_HPP_
#define _ABSTRACTMESHWRITER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>


class AbstractMeshWriter 
{
	protected:
		int mDimension; /**< Is the dimension the mesh*/
				
		std::string mPathBaseName; /**< Path to the directory where the input files are stored */
	
		std::vector< std::vector<double> > mNodeData; /**< Is an array of node coordinates ((i,j)th entry is the jth coordinate of node i)*/
		std::vector< std::vector<int> > mElementData; /**< Is an array of the nodes in each element ((i,j)th entry is the jth node of element i) */
		std::vector< std::vector<int> > mBoundaryFaceData; /**< Is an array of the nodes on each boundary face ((i,j)th entry is the jth node of face i) */		
		
		std::vector< std::vector<double> >::iterator mpNodeIterator; /**< Is an iterator for the node data */
		std::vector< std::vector<int> >::iterator mpElementIterator; /**< Is an iterator for the element data */
		std::vector< std::vector<int> >::iterator mpBoundaryFaceIterator; /**< Is an iterator for the boundary face data */		
	
		bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */
		bool mWriteMetaFile; 
	public:	
		AbstractMeshWriter() /**< Constructor */
		{
			mDimension = 0;
			mPathBaseName = "";
		}

	
	void SetNextNode(std::vector<double> nextNode);
	void SetNextElement(std::vector<int> nextElement);
	void SetNextBoundaryFace(std::vector<int> nextFace);
	void SetNextBoundaryEdge(std::vector<int> nextEdge);
	virtual void WriteFiles(){};
	int GetNumNodes(){return mNodeData.size();}
	int GetNumElements(){return mElementData.size();}
	int GetNumBoundaryFaces(){return mBoundaryFaceData.size();}
	int GetNumBoundaryEdges(){return mBoundaryFaceData.size();}
	
	virtual ~AbstractMeshWriter();
};

#endif //_ABSTRACTMESHWRITER_HPP_
