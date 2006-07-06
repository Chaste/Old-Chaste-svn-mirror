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
	if (nextNode.size() != SPACE_DIM)
	{
		throw Exception("Size of node does not match dimension.");
	}
	mNodeData.push_back(nextNode);
}

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextElement(std::vector<int> nextElement)
{
	if (nextElement.size() != ELEMENT_DIM+1)
	{
		throw Exception("Size of element does not match dimension.");
	}
	mElementData.push_back(nextElement);
}

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextBoundaryFace(std::vector<int> nextFace)
{
	if (nextFace.size() != ELEMENT_DIM)
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

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(
                                                  ConformingTetrahedralMesh<ELEMENT_DIM, 
                                                  SPACE_DIM>& rMesh)
{
    for (int i=0; i<rMesh.GetNumNodes();i++)
    {
        Node<SPACE_DIM>* p_node = rMesh.GetNodeAt(i);
        std::vector<double> coords(SPACE_DIM);
        for(int j=0; j<SPACE_DIM; j++)
        {
            coords[j] = p_node->GetPoint()[j];
        }
        SetNextNode(coords);
    }
    
    // Get an iterator over the elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter =
               rMesh.GetElementIteratorBegin();
 
    while (iter != rMesh.GetElementIteratorEnd())
    {
        std::vector<int> indices(ELEMENT_DIM+1);
        for(int j=0; j<ELEMENT_DIM+1; j++)
        {
            indices[j] = iter->GetNodeGlobalIndex(j);
        }
        SetNextElement(indices);
        iter++;
    }

    // Get a iterator over the boundary elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator boundary_iter =
               rMesh.GetBoundaryElementIteratorBegin();
    while (boundary_iter != rMesh.GetBoundaryElementIteratorEnd())
    {
        std::vector<int> indices(ELEMENT_DIM);
        for(int j=0; j<ELEMENT_DIM; j++)
        {
            indices[j] = (*boundary_iter)->GetNodeGlobalIndex(j);
        }
        SetNextBoundaryFace(indices);
        boundary_iter++;
    }
    WriteFiles();
}
 

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReader(
                                                  AbstractMeshReader<ELEMENT_DIM, 
                                                  SPACE_DIM>& rMeshReader)
{
    for (int i=0; i<rMeshReader.GetNumNodes();i++)
    {
        SetNextNode(rMeshReader.GetNextNode());
    }
    for (int i=0; i<rMeshReader.GetNumElements();i++)
    {
        SetNextElement(rMeshReader.GetNextElement());
    }
    for (int i=0; i<rMeshReader.GetNumFaces();i++)
    {
        SetNextBoundaryFace(rMeshReader.GetNextFace());
    }
    WriteFiles();
}


#endif //_ABSTRACTMESHWRITER_CPP_

