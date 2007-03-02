#ifndef _ABSTRACTMESHWRITER_CPP_
#define _ABSTRACTMESHWRITER_CPP_

#include "AbstractMeshWriter.hpp"

/**
 * Return the full path to the directory where meshes will be written.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputDirectory(void)
{
    return mpOutputFileHandler->GetTestOutputDirectory();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextNode(std::vector<double> nextNode)
{
    assert (nextNode.size() == SPACE_DIM);
    mNodeData.push_back(nextNode);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextElement(std::vector<unsigned> nextElement)
{
    assert (nextElement.size() == ELEMENT_DIM+1);
    mElementData.push_back(nextElement);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextBoundaryFace(std::vector<unsigned> nextFace)
{
    assert (nextFace.size() == ELEMENT_DIM);
    mBoundaryFaceData.push_back(nextFace);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::SetNextBoundaryEdge(std::vector<unsigned> nextEdge)
{
    SetNextBoundaryFace(nextEdge);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(
    ConformingTetrahedralMesh<ELEMENT_DIM,
    SPACE_DIM>& rMesh)
{


    NodeMap node_map(rMesh.GetNumAllNodes()); 
    unsigned new_index=0;
    for (unsigned i=0; i<(unsigned)rMesh.GetNumAllNodes();i++)
    {
        Node<SPACE_DIM>* p_node = rMesh.GetNode(i);
        
        if (p_node->IsDeleted() == false)
        {
            std::vector<double> coords(SPACE_DIM);
            for (unsigned j=0; j<SPACE_DIM; j++)
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
     assert(new_index==(unsigned)rMesh.GetNumNodes());
    
    // Get an iterator over the elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter =
        rMesh.GetElementIteratorBegin();
        
    while (iter != rMesh.GetElementIteratorEnd())
    {
        if ((*iter)->IsDeleted() == false)
        {
            std::vector<unsigned> indices(ELEMENT_DIM+1);
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                unsigned old_index=(*iter)->GetNodeGlobalIndex(j);
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
            std::vector<unsigned> indices(ELEMENT_DIM);
            for (unsigned j=0; j<ELEMENT_DIM; j++)
            {
                unsigned old_index=(*boundary_iter)->GetNodeGlobalIndex(j);
                indices[j] = node_map.GetNewIndex(old_index); 
            }
            SetNextBoundaryFace(indices);
        }
        boundary_iter++;
    }
    WriteFiles();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReader(
    AbstractMeshReader<ELEMENT_DIM,
    SPACE_DIM>& rMeshReader)
{
    for (unsigned i=0; i<rMeshReader.GetNumNodes();i++)
    {
        SetNextNode(rMeshReader.GetNextNode());
    }
    for (unsigned i=0; i<rMeshReader.GetNumElements();i++)
    {
        SetNextElement(rMeshReader.GetNextElement());
    }
    for (unsigned i=0; i<rMeshReader.GetNumFaces();i++)
    {
        SetNextBoundaryFace(rMeshReader.GetNextFace());
    }
    WriteFiles();
}


#endif //_ABSTRACTMESHWRITER_CPP_

