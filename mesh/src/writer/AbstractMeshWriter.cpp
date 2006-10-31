#ifndef _ABSTRACTMESHWRITER_CPP_
#define _ABSTRACTMESHWRITER_CPP_

#include "AbstractMeshWriter.hpp"

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
    assert (nextNode.size() == SPACE_DIM);
    mNodeData.push_back(nextNode);
}

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextElement(std::vector<int> nextElement)
{
    assert (nextElement.size() == ELEMENT_DIM+1);
    mElementData.push_back(nextElement);
}

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextBoundaryFace(std::vector<int> nextFace)
{
    assert (nextFace.size() == ELEMENT_DIM);
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


    NodeMap node_map(rMesh.GetNumAllNodes()); 
    int new_index=0;
    for (int i=0; i<rMesh.GetNumAllNodes();i++)
    {
        Node<SPACE_DIM>* p_node = rMesh.GetNodeAt(i);
        
        if (p_node->IsDeleted() == false)
        {
            std::vector<double> coords(SPACE_DIM);
            for (int j=0; j<SPACE_DIM; j++)
            {
                coords[j] = p_node->GetPoint()[j];
            }
            SetNextNode(coords);
            node_map.SetNewIndex(i,new_index++);
        }
        else
        {
            node_map.SetDeleted(i);
        }
    }
    assert(new_index==rMesh.GetNumNodes());
    
    // Get an iterator over the elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter =
        rMesh.GetElementIteratorBegin();
        
    while (iter != rMesh.GetElementIteratorEnd())
    {
        if ((*iter)->IsDeleted() == false)
        {
            std::vector<int> indices(ELEMENT_DIM+1);
            for (int j=0; j<ELEMENT_DIM+1; j++)
            {
                int old_index=(*iter)->GetNodeGlobalIndex(j);
                indices[j] = node_map.GetNewIndex(old_index); 
            }
            SetNextElement(indices);
        }
        
        iter++;
    }
    
    // Get a iterator over the boundary elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator boundary_iter =
        rMesh.GetBoundaryElementIteratorBegin();
    while (boundary_iter != rMesh.GetBoundaryElementIteratorEnd())
    {
        if ((*boundary_iter)->IsDeleted() == false)
        {
            std::vector<int> indices(ELEMENT_DIM);
            for (int j=0; j<ELEMENT_DIM; j++)
            {
                int old_index=(*boundary_iter)->GetNodeGlobalIndex(j);
                indices[j] = node_map.GetNewIndex(old_index); 
            }
            SetNextBoundaryFace(indices);
        }
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

