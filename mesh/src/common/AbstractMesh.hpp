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
#include "BoundaryElement.hpp"
#include "Element.hpp"
#include "AbstractMeshReader.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractMesh
{
protected:  // Give access of these variables to subclasses
    std::vector<Node<SPACE_DIM> *> mNodes;
    std::vector<Node<SPACE_DIM> *> mBoundaryNodes;

    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> mElements;
    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;
    
    std::vector<unsigned> mNodesPerProcessor;   

public:
    typedef typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::const_iterator ElementIterator;
    typedef typename std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *>::const_iterator BoundaryElementIterator;
    typedef typename std::vector<Node<SPACE_DIM> *>::const_iterator BoundaryNodeIterator;

    virtual ~AbstractMesh()
    {
    }

    virtual unsigned GetNumNodes() const;
    virtual unsigned GetNumElements() const;
    virtual unsigned GetNumBoundaryElements() const;
    unsigned GetNumBoundaryNodes();// should this be overloaded and virtual too?

    unsigned GetNumAllNodes() const;
    unsigned GetNumAllElements();
    unsigned GetNumAllBoundaryElements();

    Node<SPACE_DIM> *GetNode(unsigned index);    
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index);
    
    /**
     * Sets the ownership of each element according to which nodes are owned by the
     * process.
     * @param lo is the lowest node number owned by the process
     * @param hi is one higher than the highest node number owned by the process
     * ie. this process owns nodes [lo..hi)
     * and element is "owned" if one or more of its nodes are owned
     */
    virtual void SetElementOwnerships(unsigned lo, unsigned hi)=0;
    
    virtual void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                         bool cullInternalFaces=false)=0;
    
    virtual void ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile);

    std::vector<unsigned>& rGetNodesPerProcessor();
    
    virtual void PermuteNodes();
    
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

private:
    virtual unsigned SolveNodeMapping(unsigned index)=0;
    virtual unsigned SolveElementMapping(unsigned index)=0;        

};

/// Returns the number of nodes that are actually in use
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return this->mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryNodes()
{
    return this->mBoundaryNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes() const
{
    return this->mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements()
{
    return this->mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllBoundaryElements()
{
    return this->mBoundaryElements.size();
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return this->mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index)
{
    unsigned local_index = SolveNodeMapping(index);
    return this->mNodes[local_index];
}
    
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index)
{
    unsigned local_index = SolveElementMapping(index);
    return this->mElements[local_index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile)
{
    NEVER_REACHED
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned>& AbstractMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodesPerProcessor()
{
    return mNodesPerProcessor;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    NEVER_REACHED
}

#endif /*ABSTRACTMESH_HPP_*/
