#ifndef _ABSTRACTMESHWRITER_CPP_
#define _ABSTRACTMESHWRITER_CPP_

#include "AbstractMeshWriter.hpp"
#include "Exception.hpp"

/**
 * Return the full path to the directory where meshes will be written.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
std::string AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputDirectory(void)
{
    return mpOutputFileHandler->GetTestOutputDirectory();
}

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextNode(std::vector<double> nextNode)
{
	if (nextNode.size() != mDimension)
	{
		throw Exception("Size of node does not match dimension.");
	}
	mNodeData.push_back(nextNode);
}

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextElement(std::vector<int> nextElement)
{
	if (nextElement.size() != mDimension+1)
	{
		throw Exception("Size of element does not match dimension.");
	}
	mElementData.push_back(nextElement);
}

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextBoundaryFace(std::vector<int> nextFace)
{
	if (nextFace.size() != mDimension)
	{
		throw Exception("Size of face or edge does not match dimension.");
	}
	mBoundaryFaceData.push_back(nextFace);
}

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextBoundaryEdge(std::vector<int> nextEdge)
{
	SetNextBoundaryFace(nextEdge);
}

#endif //_ABSTRACTMESHWRITER_CPP_

