#ifndef _ABSTRACTMESHWRITER_HPP_
#define _ABSTRACTMESHWRITER_HPP_

#include "AbstractMeshReadWrite.hpp"

class AbstractMeshWriter : public AbstractMeshReadWrite
{
public:
	
	void SetNextNode(std::vector<double> nextNode);
	void SetNextElement(std::vector<int> nextElement);
	void SetNextBoundaryFace(std::vector<int> nextFace);
	void SetNextBoundaryEdge(std::vector<int> nextEdge);
	
	int GetNumNodes(){return mNodeData.size();}
	int GetNumElements(){return mElementData.size();}
	int GetNumBoundaryFaces(){return mBoundaryFaceData.size();}
	int GetNumBoundaryEdges(){return mBoundaryFaceData.size();}
	
	virtual ~AbstractMeshWriter();
};

#endif //_ABSTRACTMESHWRITER_HPP_
