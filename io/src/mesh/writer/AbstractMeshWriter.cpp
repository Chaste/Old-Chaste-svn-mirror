#include "AbstractMeshWriter.hpp"
#include "global/src/Exception.hpp"


void AbstractMeshWriter::SetNextNode(std::vector<double> nextNode)
{
	if (nextNode.size() != mDimension)
	{
		throw Exception("Size of node does not match dimension.");
	}
	mNodeData.push_back(nextNode);
}

void AbstractMeshWriter::SetNextElement(std::vector<int> nextElement)
{
	if (nextElement.size() != mDimension+1)
	{
		throw Exception("Size of element does not match dimension.");
	}
	mElementData.push_back(nextElement);
}

void AbstractMeshWriter::SetNextBoundaryFace(std::vector<int> nextFace)
{
	if (nextFace.size() != mDimension)
	{
		throw Exception("Size of face or edge does not match dimension.");
	}
	mBoundaryFaceData.push_back(nextFace);
}

void AbstractMeshWriter::SetNextBoundaryEdge(std::vector<int> nextEdge)
{
	SetNextBoundaryFace(nextEdge);
}

AbstractMeshWriter::~AbstractMeshWriter()
{
}
