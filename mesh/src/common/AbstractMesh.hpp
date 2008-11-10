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
#ifndef ABSTRACTMESH_HPP_
#define ABSTRACTMESH_HPP_

#include "Node.hpp"
#include "Element.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM,
         class NODE_CONTAINER, class ELEMENT_CONTAINER, class BOUNDARY_CONTAINER>
class AbstractMesh
{
public:
    typedef typename ELEMENT_CONTAINER::const_iterator ElementIterator;
    typedef typename BOUNDARY_CONTAINER::const_iterator BoundaryElementIterator;
    typedef typename NODE_CONTAINER::const_iterator BoundaryNodeIterator;

protected:  // Give access of these variables to subclasses
    NODE_CONTAINER mNodes;
    NODE_CONTAINER mBoundaryNodes;

    ELEMENT_CONTAINER mElements;
    BOUNDARY_CONTAINER mBoundaryElements;   
    
public:

    virtual unsigned GetNumNodes();
    virtual unsigned GetNumElements();
    virtual unsigned GetNumBoundaryElements();
    unsigned GetNumBoundaryNodes();// should this be overloaded and virtual too?

    virtual bool GetNodeIsLocal(unsigned index)=0;
    virtual bool GetElementIsLocal(unsigned index)=0;    
    
    Node<SPACE_DIM> *GetNode(unsigned index);    
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index);
    
    
/// \todo: move implementations out of class definition    
    
    /**
     * Return a pointer to the first element in the mesh.
     */
    ElementIterator GetElementIteratorBegin() const
    {
        return mElements.begin();
    }
    /**
     * Return a pointer to *one past* the last element in the mesh
     * (for consistency with STL iterators).
     */
    ElementIterator GetElementIteratorEnd() const
    {
        return mElements.end();
    }

    /**
     * Return a pointer to the first boundary element in the mesh.
     */

    BoundaryElementIterator GetBoundaryElementIteratorBegin() const
    {
        return mBoundaryElements.begin();
    }
    /**
     * Return a pointer to *one past* the last boundary element in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryElementIterator GetBoundaryElementIteratorEnd() const
    {
        return mBoundaryElements.end();
    }

    /**
     * Return a pointer to the first boundary node in the mesh.
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorBegin() const
    {
        return mBoundaryNodes.begin();
    }
    /**
     * Return a pointer to *one past* the last boundary node in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorEnd() const
    {
        return mBoundaryNodes.end();
    }
    
};

/// Returns the number of nodes that are actually in use
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class NODE_CONTAINER, class ELEMENT_CONTAINER, class BOUNDARY_CONTAINER>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM, NODE_CONTAINER, ELEMENT_CONTAINER, BOUNDARY_CONTAINER>::GetNumNodes()
{
    return this->mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class NODE_CONTAINER, class ELEMENT_CONTAINER, class BOUNDARY_CONTAINER>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM, NODE_CONTAINER, ELEMENT_CONTAINER, BOUNDARY_CONTAINER>::GetNumElements()
{
    return this->mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class NODE_CONTAINER, class ELEMENT_CONTAINER, class BOUNDARY_CONTAINER>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM, NODE_CONTAINER, ELEMENT_CONTAINER, BOUNDARY_CONTAINER>::GetNumBoundaryNodes()
{
    return this->mBoundaryNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class NODE_CONTAINER, class ELEMENT_CONTAINER, class BOUNDARY_CONTAINER>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM, NODE_CONTAINER, ELEMENT_CONTAINER, BOUNDARY_CONTAINER>::GetNumBoundaryElements()
{
    return this->mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class NODE_CONTAINER, class ELEMENT_CONTAINER, class BOUNDARY_CONTAINER>
Node<SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM, NODE_CONTAINER, ELEMENT_CONTAINER, BOUNDARY_CONTAINER>::GetNode(unsigned index)
{
    assert(GetNodeIsLocal(index));
    return this->mNodes[index];
}
    
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class NODE_CONTAINER, class ELEMENT_CONTAINER, class BOUNDARY_CONTAINER>
Element<ELEMENT_DIM, SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM, NODE_CONTAINER, ELEMENT_CONTAINER, BOUNDARY_CONTAINER>::GetElement(unsigned index)
{
    assert(GetElementIsLocal(index));
    return this->mElements[index];
}


#endif /*ABSTRACTMESH_HPP_*/
