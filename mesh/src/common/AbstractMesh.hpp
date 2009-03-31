/*

Copyright (C) University of Oxford, 2005-2009

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

/** 
 * Abstract base class for all meshes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractMesh
{
private:

    /**
     * Pure virtual solve node mapping method.
     * Overridden in TetrahedralMesh and ParallelTetrahedralMesh classes.
     * 
     * @param index
     */
    virtual unsigned SolveNodeMapping(unsigned index) const = 0;

    /**
     * Pure virtual solve element mapping method.
     * Overridden in TetrahedralMesh and ParallelTetrahedralMesh classes.
     * 
     * @param index
     */
    virtual unsigned SolveElementMapping(unsigned index) const = 0;

    /**
     * Pure virtual solve boundary element mapping method.
     * Overridden in TetrahedralMesh and ParallelTetrahedralMesh classes.
     * 
     * @param index
     */
    virtual unsigned SolveBoundaryElementMapping(unsigned index) const = 0;

protected:  // Give access of these variables to subclasses

    /** Vector of pointers to nodes in the mesh. */
    std::vector<Node<SPACE_DIM> *> mNodes;

    /** Vector of pointers to boundary nodes in the mesh. */
    std::vector<Node<SPACE_DIM> *> mBoundaryNodes;

    /** Vector of pointers to elements in the mesh. */
    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> mElements;

    /** Vector of pointers to boundary elements in the mesh. */
    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;

    /** Vector containing the number of nodes owned by each processor. */
    std::vector<unsigned> mNodesPerProcessor;

    /** Vector containing node permutation information. */
    std::vector<unsigned> mNodesPermutation;

    /** 
     * If the mesh is constructed from file using a MeshReader, this member  
     * variable stores the base name of these files. 
     */
    std::string mMeshFileBaseName;

public:

    /** Definition of element Iterator type. */
    typedef typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::const_iterator ElementIterator;
    /** Definition of boundary element Iterator type. */
    typedef typename std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *>::const_iterator BoundaryElementIterator;
    /** Definition of boundary node Iterator type. */
    typedef typename std::vector<Node<SPACE_DIM> *>::const_iterator BoundaryNodeIterator;

    /** 
     * Constructor.
     */
    AbstractMesh();

    /** 
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractMesh();

    /** 
     * Get the number of nodes that are actually in use.
     */
    virtual unsigned GetNumNodes() const;

    /** 
     * Get the number of elements that are actually in use.
     */
    virtual unsigned GetNumElements() const;

    /** 
     * Get the number of boundary elements that are actually in use.
     */
    virtual unsigned GetNumBoundaryElements() const;

    /** 
     * Get the number of boundary nodes in the mesh.
     */
    unsigned GetNumBoundaryNodes();// should this be overloaded and virtual too?

    /** 
     * Get the total number of nodes (including those marked as deleted).
     */
    unsigned GetNumAllNodes() const;

    /** 
     * Get the total number of elements (including those marked as deleted).
     */
    unsigned GetNumAllElements();

    /** 
     * Get the total number of boundary elements (including those marked as deleted).
     */
    unsigned GetNumAllBoundaryElements();

    /**  
     * Get the node with a given index in the mesh.
     * 
     * @param index
     * @return a pointer to the node.
     */
    Node<SPACE_DIM> *GetNode(unsigned index) const;

    /**  
     * Get the element with a given index in the mesh.
     * 
     * @param index
     * @return a pointer to the element.
     */
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**  
     * Get the boundary element with a given index in the mesh.
     * 
     * @param index
     * @return a pointer to the boundary element.
     */
    BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* GetBoundaryElement(unsigned index) const;

    /**
     * Sets the ownership of each element according to which nodes are owned by the
     * process.
     * 
     * @param lo is the lowest node number owned by the process
     * @param hi is one higher than the highest node number owned by the process
     * ie. this process owns nodes [lo..hi)
     * and element is "owned" if one or more of its nodes are owned
     */
    virtual void SetElementOwnerships(unsigned lo, unsigned hi);

    /** 
     * Construct the mesh using a MeshReader. 
     * This method must be overridden in concrete classes. 
     *  
     * @param rMeshReader the mesh reader 
     * @param cullInternalFaces whether to cull internal faces (defaults to false) 
     */ 
    virtual void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                         bool cullInternalFaces=false)=0;

    /** 
     * Read in the number of nodes per processor from file. 
     *  
     * @param nodesPerProcessorFile 
     */ 
    virtual void ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile);

    /** 
     * Get method for mNodesPerProcessor.
     */
    std::vector<unsigned>& rGetNodesPerProcessor();

    /** 
     * Permute the nodes so that they appear in a different order in mNodes 
     * (and their mIndex's are altered accordingly).
     */
    virtual void PermuteNodes();      

    /**
     * Return a pointer to the first element in the mesh.
     */
    ElementIterator GetElementIteratorBegin() const;

    /**
     * Return a pointer to *one past* the last element in the mesh
     * (for consistency with STL iterators).
     */
    ElementIterator GetElementIteratorEnd() const;

    /**
     * Return a pointer to the first boundary element in the mesh.
     */
    BoundaryElementIterator GetBoundaryElementIteratorBegin() const;

    /**
     * Return a pointer to *one past* the last boundary element in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryElementIterator GetBoundaryElementIteratorEnd() const;

    /**
     * Return a pointer to the first boundary node in the mesh.
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorBegin() const;

    /**
     * Return a pointer to *one past* the last boundary node in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorEnd() const;

    /**
     * Compute the inverse Jacobian for a given element in the mesh.
     * 
     * @param elementIndex index of an element
     * @param rJacobian  the Jacobian matrix
     * @param rJacobianDeterminant  the determinant of the Jacobian matrix
     * @param rInverseJacobian  the inverse Jacobian matrix
     */
    virtual void GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double &rJacobianDeterminant, c_matrix<double, SPACE_DIM, SPACE_DIM>& rInverseJacobian) const;

    /**
     * Compute the weighted direction for a given boundary element.
     * 
     * @param elementIndex index of an element
     * @param rWeightedDirection the weighted direction vector
     * @param rJacobianDeterminant  the determinant of the Jacobian matrix
     */
    virtual void GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant) const;

    /** 
     * Get method for mMeshFileBaseName.
     */
    std::string GetMeshFileBaseName() const;

    /** 
     * Get method for mNodesPermutation.
     */
    std::vector<unsigned>& rGetNodePermutation();
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
{
    assert(hi>=lo);
    for (unsigned element_index=0; element_index<this->mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element=this->mElements[element_index];
        p_element->SetOwnership(false);
        for (unsigned local_node_index=0; local_node_index< p_element->GetNumNodes(); local_node_index++)
        {
            unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
            if (lo<=global_node_index && global_node_index<hi)
            {
                p_element->SetOwnership(true);
                break;
            }
        }

    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMesh<ELEMENT_DIM, SPACE_DIM>::AbstractMesh()
: mMeshFileBaseName("")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMesh<ELEMENT_DIM, SPACE_DIM>::~AbstractMesh()
{
    // Iterate over nodes and free the memory
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
    // Iterate over elements and free the memory
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        delete this->mElements[i];
    }
    // Iterate over boundary elements and free the memory
    for (unsigned i=0; i<this->mBoundaryElements.size(); i++)
    {
        delete this->mBoundaryElements[i];
    }
}

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
Node<SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index) const
{
    unsigned local_index = SolveNodeMapping(index);
    return this->mNodes[local_index];
}
    
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    unsigned local_index = SolveElementMapping(index);
    return this->mElements[local_index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElement(unsigned index) const
{
    unsigned local_index = SolveBoundaryElementMapping(index);
    return this->mBoundaryElements[local_index];
}    

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile)
{
    NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned>& AbstractMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodesPerProcessor()
{
    return mNodesPerProcessor;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    NEVER_REACHED;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin() const
{
    return mElements.begin();
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd() const
{
    return mElements.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorBegin() const
{
    return mBoundaryElements.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorEnd() const
{
    return mBoundaryElements.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryNodeIteratorBegin() const
{
    return mBoundaryNodes.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryNodeIteratorEnd() const
{
    return mBoundaryNodes.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double &rJacobianDeterminant, c_matrix<double, SPACE_DIM, SPACE_DIM>& rInverseJacobian) const
{
    mElements[SolveElementMapping(elementIndex)]->CalculateInverseJacobian(rJacobian, rJacobianDeterminant, rInverseJacobian);
}    

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant) const
{
    mBoundaryElements[SolveBoundaryElementMapping(elementIndex)]->CalculateWeightedDirection(rWeightedDirection, rJacobianDeterminant );
}    

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetMeshFileBaseName() const
{
    if (mMeshFileBaseName == "")
    {
        EXCEPTION("This mesh was not constructed from a file.");
    }
    
    return mMeshFileBaseName;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned>& AbstractMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodePermutation()
{
    return mNodesPermutation;
}

#endif /*ABSTRACTMESH_HPP_*/
