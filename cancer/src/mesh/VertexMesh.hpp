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


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh : public AbstractMesh< ELEMENT_DIM, SPACE_DIM>
{
private:
    
    unsigned SolveNodeMapping(unsigned index) const;
    unsigned SolveElementMapping(unsigned index) const;
    unsigned SolveBoundaryElementMapping(unsigned index) const;   
    
    std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> mVertexElements;
    
public:
    
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                         bool cullInternalFaces=false)
    {};
    
    void SetElementOwnerships(unsigned lo, unsigned hi)
    {};
    
    VertexMesh(std::vector<Node<SPACE_DIM> *> nodes);
    
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertex_elements);
    
    unsigned GetNumVertexElements();
    
    void Clear();
    
    ~VertexMesh();
    
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM> *> nodes)
{
    Clear();
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<SPACE_DIM>* temp_node = nodes[index];
        this->mNodes.push_back(temp_node);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, 
                        std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertex_elements)
{
    Clear();
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<SPACE_DIM>* temp_node = nodes[index];
        this->mNodes.push_back(temp_node);
    }
    
    for (unsigned index=0; index<vertex_elements.size(); index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* temp_vertex_element = vertex_elements[index];
        this->mVertexElements.push_back(temp_vertex_element);
    }
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    for (unsigned i=0; i<this->mBoundaryElements.size(); i++)
    {
        delete this->mBoundaryElements[i];
    }
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        delete this->mElements[i];
    }
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }

    this->mNodes.clear();
    this->mElements.clear();
    this->mBoundaryElements.clear();
    this->mBoundaryNodes.clear();      
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    assert(index < this->mBoundaryElements.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumVertexElements()
{
    return mVertexElements.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::~VertexMesh()
{};
    

#endif /*VERTEXMESH_HPP_*/
