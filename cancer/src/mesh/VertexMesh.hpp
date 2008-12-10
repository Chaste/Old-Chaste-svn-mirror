/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef VERTEXMESH_HPP_
#define VERTEXMESH_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>


#include "AbstractMesh.hpp"
#include "VertexElement.hpp"
#include "BoundaryElement.hpp"
#include "NodeMap.hpp"
#include "Node.hpp"
#include "Exception.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MutableMesh.hpp"
#include "VoronoiTessellation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh
{
private:
    std::vector<Node<SPACE_DIM> *> mNodes;
    std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> mElements;
    
    void SetupVertexElementsOwnedByNodes();
    
public:
    /**
     *  Constructor takes in node and VertexElements
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertex_elements);
    ~VertexMesh();
    
    VertexMesh(unsigned numAcross,unsigned numUp);

    unsigned GetNumNodes() const;
    unsigned GetNumElements() const;

//// when will these be needed?
//    unsigned GetNumAllNodes() const;
//    unsigned GetNumAllElements();

    Node<SPACE_DIM>* GetNode(unsigned index) const;    
    VertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;


//    
//    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
//                                         bool cullInternalFaces=false)
//    {}
    
    std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> GetElementsOwnedByNode(Node<SPACE_DIM>* p_node);
    
    void Clear();
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, 
                                               std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertex_elements)
{
    Clear();
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<SPACE_DIM>* temp_node = nodes[index];
        mNodes.push_back(temp_node);
    }
    
    for (unsigned index=0; index<vertex_elements.size(); index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* temp_vertex_element = vertex_elements[index];
        mElements.push_back(temp_vertex_element);
    }
    
    SetupVertexElementsOwnedByNodes();
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(unsigned numAcross,unsigned numUp)
{
    HoneycombMeshGenerator generator(numAcross+1,numUp+1,0,false);
    MutableMesh<2,2>* p_mesh = generator.GetMesh();
    VoronoiTessellation<2> tessellation(*p_mesh);
    
    for (unsigned i = 0;i<tessellation.GetNumVertices();i++)
    {
        c_vector<double,2>* position = tessellation.GetVertex(i);
        Node<2>* p_node = new Node<2>(0, false, (*position)(0), (*position)(1));
        mNodes.push_back(p_node);
    }    

    // todo: loop over the p_mesh's nodes, and if it is a non-boundary node create a VertexElement using
    // the corresponding cell. Then get rid of the nodes in mNodes that do not belong in any cell.

}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetupVertexElementsOwnedByNodes()
{
    for (unsigned index=0; index<mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_temp_vertex_element = mElements[index];
        for (unsigned node_index=0; node_index<p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        delete mNodes[i];
    }

    mNodes.clear();
    mElements.clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNodes.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index) const
{
    assert(index < mNodes.size());
    return mNodes[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM,SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}


//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes()
//{
//    return mNodes.size();
//}
//
//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements()
//{
//    return mElements.size();
//}
//


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::~VertexMesh()
{
    //Clear();
}
    

#endif /*VERTEXMESH_HPP_*/
