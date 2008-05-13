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


#ifndef _CONFORMINGTETRAHEDRALMESH_HPP_
#define _CONFORMINGTETRAHEDRALMESH_HPP_

#include <boost/serialization/access.hpp>

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>


//Jonathan Shewchuk's triangle
#define REAL double
#define VOID void
#include "triangle.h"
#undef REAL

#include "AbstractMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "Element.hpp"
#include "BoundaryElement.hpp"
#include "Node.hpp"
#include "NodeMap.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"

#ifndef SPECIAL_SERIAL
#include "PetscTools.hpp"
#endif //SPECIAL_SERIAL

#include <boost/serialization/export.hpp>

//////////////////////////////////////////////////////////////////////////
//   DECLARATION
//////////////////////////////////////////////////////////////////////////
 
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ConformingTetrahedralMesh
{
    friend class TestConformingTetrahedralMesh; // to give access to private methods (not variables)
    friend class TestCryptSimulation2d; // to give access to private methods (not variables)
public:
    typedef typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::const_iterator ElementIterator;
    typedef typename std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *>::const_iterator BoundaryElementIterator;
    typedef typename std::vector<Node<SPACE_DIM> *>::const_iterator BoundaryNodeIterator;

protected:  // Give access of these variables to subclasses
    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> mElements;
    std::vector<Node<SPACE_DIM> *> mNodes;
    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;
    
    std::vector<unsigned> mNodesPerProcessor;
    
    /// Indices of elements/nodes that have been deleted - these indices can be reused when adding
    /// new elements/nodes
    std::vector<unsigned> mDeletedElementIndices;
    std::vector<unsigned> mDeletedBoundaryElementIndices;
    std::vector<unsigned> mDeletedNodeIndices;
    bool mAddedNodes;
    std::vector< Node<SPACE_DIM> *> mBoundaryNodes;
    
private:    
    //unsigned mNumCornerNodes;
   
    std::map<unsigned, unsigned> mSmasrmIndexMap;
    
    /**
     * Check whether any neighbouring node is inside the circumsphere of this element.
     * @param pointer to an element
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of the element, as a proportion of the circumsphere radius.
     */
    bool CheckVoronoi(Element<ELEMENT_DIM, SPACE_DIM>  *pElement, double maxPenetration);
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       // Don't do anything - this is just so subclasses can archive member variables.
    }
    
public:

    ConformingTetrahedralMesh();
    ConformingTetrahedralMesh(unsigned numElements);
    ConformingTetrahedralMesh(std::vector<Node<SPACE_DIM> *> nodes);
    
    virtual ~ConformingTetrahedralMesh();
    
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                 bool cullInternalFaces=false);
    
    void RescaleMeshFromBoundaryNode(ChastePoint<1> updatedPoint, unsigned boundaryNodeIndex);
    
    Node<SPACE_DIM> *GetNode(unsigned index);
    
    unsigned GetNumNodes();
    unsigned GetNumElements();
    unsigned GetNumBoundaryElements();
    unsigned GetNumAllNodes() const;
    unsigned GetNumAllElements();
    unsigned GetNumAllBoundaryElements();
    unsigned GetNumBoundaryNodes();
    
    void ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile);
    std::vector<unsigned>& rGetNodesPerProcessor();
    
    /** 
     * Add a node to the mesh.
     * 
     * NB. After calling this one or more times, you must then call ReMesh
     *
     */
    
    virtual unsigned AddNode(Node<SPACE_DIM> *pNewNode);
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
    
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index)
    {
        return (mElements[index]);
    }
    
    BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* GetBoundaryElement(unsigned index)
    {
        return (mBoundaryElements[index]);
    }
    
    virtual void SetNode(unsigned index, ChastePoint<SPACE_DIM> point, bool concreteMove=true);
    void MoveMergeNode(unsigned index, unsigned targetIndex, bool concreteMove=true);

    void DeleteNode(unsigned index);
    
    void DeleteNodePriorToReMesh(unsigned index);
    
    unsigned RefineElement(Element<ELEMENT_DIM,SPACE_DIM>* pElement, ChastePoint<SPACE_DIM> Point);
    void RefreshMesh(void);
    
    /**
     * Return the volume of a mesh, calculated by adding the determinant of each element 
     * and dividing by n!, where n is the element dimension.
     */
    double CalculateVolume();
    double CalculateSurfaceArea();
    
    void Translate(c_vector<double, SPACE_DIM> displacement);
    void Translate(const double xMovement=0.0, const double yMovement=0.0, const double zMovement=0.0);
    void Scale(const double xFactor=1.0, const double yFactor=1.0, const double zFactor=1.0);
    void Rotate(c_matrix<double , SPACE_DIM, SPACE_DIM> rotation_matrix);
    void Rotate(c_vector<double,3> axis, double angle);
    void RotateX(const double theta);
    void RotateY(const double theta);
    void RotateZ(const double theta);
    /**Rotating a 2D mesh equates that rotation around the z-axis*/
    void Rotate(double theta)
    {
        RotateZ(theta);
    }
    
    /**
     * Remove a boundary node, and update all the appropriate data structures.
     * 
     * The deleted node is not removed from the list, merely marked as deleted,
     * and can be reused when a new node is added to the mesh.
     * 
     * Any elements or boundary elements containing this node will be removed.
     * The boundary nodes information will be updated with new boundary node(s).
     * NB: New boundary elements WILL NOT be added.
     * 
     * @param index  The index of the node to remove.
     */
    void DeleteBoundaryNodeAt(unsigned index);
    
    /**
     * Re-index a mesh so that it has no deleted elements or nodes
     */
    void ReIndex(NodeMap& map);
    
    /**
     * Re-mesh a mesh using triangle or tetgen
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     */
    virtual void ReMesh(NodeMap& map);

    /** 
     * Alternative version of remesh which takes no parameters does not require a NodeMap. Note: inherited
     * classes should overload ReMesh(NodeMap&)
     */
    void ReMesh();
    
    /**
     * Re-mesh a mesh using triangle in 2D via library calls
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     */
    virtual void ReMeshWithTriangleLibrary(NodeMap& map);
    
    /**
     * Permute the nodes so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly).
     * 
     */
    void PermuteNodes();
    
    /**
      * Permute the nodes so that they appear in a different order in mNodes
      * (and their mIndex's are altered accordingly) using Metis binaries.
      * 
      * @param numProcs Number of processors (e.g. number of partitions)
      */
    void PermuteNodesWithMetisBinaries(unsigned numProcs);
    
    /**
      * Permute the nodes so that they appear in a different order in mNodes
      * (and their mIndex's are altered accordingly).
     * @param perm is a vector containing the new indices
     */
    void PermuteNodes(std::vector<unsigned>& perm);
    
    /**
     * Checks the entire mesh element by element and checks whether any neighbouring node
     * is inside the circumsphere of this element.
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of an element that it is not a member of, as a proportion of the
     * circumsphere radius.
     */
    bool CheckVoronoi(double maxPenetration=0.0);
    
    void ConstructLinearMesh(unsigned width);
    
    /**
     * Construct a rectangular grid on [0,width]x[0,height]
     * diagonals can be staggered so that there is no prefered diffusion propagation
     * direction.
     */
    void ConstructRectangularMesh(unsigned width, unsigned height, bool stagger=true);
    
    /**
     * Construct a cuboid grid on [0,width]x[0,height]x[0,depth]
     * diagonals can be staggered so that there is no prefered diffusion propagation
     * direction.
     * 
     * @param stagger whether the mesh should 'jumble' up the elements (defaults to false)
     */
    void ConstructCuboid(unsigned width, unsigned height, unsigned depth, bool stagger=false);
    
    /**
     *  Returns the element index for the first element that is known to contain a test point
     *  @param testPoint
     *  @param strict Should the element returned contain the point in the interior and
     *  not on an edge/face/vertex (default = not strict)
     *  @param a set of guesses for the element (a set of element indices), to be checked
     *  first for potential efficiency improvements. (default = empty set)
     */
    unsigned GetContainingElementIndex(ChastePoint<SPACE_DIM> testPoint, bool strict=false, std::set<unsigned> testElements=std::set<unsigned>());
    
    /**
     *  Returns the element index for an element is closest to the testPoint
     * "Closest" means that the minimum interpolation weights for the testPoint are 
     * maximised for this element
     *  @param testPoint
     * 
     */
    unsigned GetNearestElementIndex(ChastePoint<SPACE_DIM> testPoint);
    
    /**
     *  Returns all element indices for elements that are known to contain a test point
     *  @param testPoint 
     */
    std::vector<unsigned> GetContainingElementIndices(ChastePoint<SPACE_DIM> testPoint);
    
    
    /**
     * Sets the ownership of each element according to which nodes are owned by the
     * process.
     * @param lo is the lowest node number owned by the process
     * @param hi is one higher than the highest node number owned by the process
     * ie. this process owns nodes [lo..hi)
     * and element is "owned" if one or more of its nodes are owned
     */
    void SetElementOwnerships(unsigned lo, unsigned hi);
    
    /**
     *  Clear all the data in the mesh
     */
    void Clear();
    
    /**
     *  Return the set of nodes which are on the boundary of the flagged region(s)
     */  
    std::set<unsigned> CalculateBoundaryOfFlaggedRegion();
    
    /**
     * Returns a vector between two points in space
     * 
     * @param rLocationA a c_vector of co-ordinates
     * @param rLocationB a c_vector of co-ordinates
     * 
     * @return vector from location A to location B.
     * 
     * N.B. This can be overwritten in daughter classes
     * e.g. Cylindrical2dMesh.
     * 
     */
    virtual c_vector<double, SPACE_DIM> GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB);
    
    /**
     * Calcuates the angle between the node at indexB and the x axis about
     * the node at indexA. The angle returned is in the range (-pi,pi]
     */
    double GetAngleBetweenNodes(unsigned indexA, unsigned indexB);
    
    /**
     * Calculates the `width' of any dimension of the mesh.
     * 
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     * 
     * N.B. Overwritten in Cylindrical2dMesh
     */    
    virtual double GetWidth(const unsigned& rDimension) const;
    
    /**
     * Calculates the `width extremes' of any dimension of the mesh.
     * 
     * @param rDimension a dimension (0,1 or 2)
     * @return The minimum and maximum co-ordinates of any node in this dimension.
     * 
     */    
    c_vector<double,2> GetWidthExtremes(const unsigned& rDimension) const;
    
    void UnflagAllElements();
    
    
    /**
     *  Flag all elements not containing ANY of the given nodes
     */
    void FlagElementsNotContainingNodes(std::set<unsigned> nodesList);
    
    /**
     *  Set up a map between the nodes in the flagged region of the mesh and 
     */
    void SetupSmasrmMap()
    {
        // Figure out the SMASRM size, and generate a map from global node number
        // to SMASRM index.
        mSmasrmIndexMap.clear();

        ElementIterator iter = GetElementIteratorBegin();
        unsigned smasrm_size = 0;

        while (iter != GetElementIteratorEnd())
        {
            Element<ELEMENT_DIM, SPACE_DIM>& element = **iter;
            
            if (element.IsFlagged())
            {
                // Add this element's nodes to the map
                const unsigned num_nodes = element.GetNumNodes();
                for (unsigned i=0; i<num_nodes; i++)
                {
                    unsigned node_index = element.GetNodeGlobalIndex(i);
                    if (mSmasrmIndexMap.count(node_index) == 0)
                    {
                        // This is a new node
                        mSmasrmIndexMap[node_index] = smasrm_size++;
                    }
                }
            }
            ++iter;
        }
        assert(mSmasrmIndexMap.size() == smasrm_size);
    }
     
     
    /** 
     *  Get the map between the flagged node indices and
     */
    std::map<unsigned, unsigned>& rGetSmasrmMap()
    {
        return mSmasrmIndexMap;
    }

    /**
     * Iterator over edges in the mesh, which correspond to springs between cells.
     * 
     * This class takes care of the logic to make sure that you consider each edge exactly once.
     */
    class EdgeIterator
    {
    public:
        /**
         * Get a pointer to the node in the mesh at end A of the spring.
         */
        Node<SPACE_DIM>* GetNodeA();
        /**
         * Get a pointer to the node in the mesh at end B of the spring.
         */
        Node<SPACE_DIM>* GetNodeB();
                
        bool operator!=(const EdgeIterator& other);
        
        /**
         * Prefix increment operator.
         */
        EdgeIterator& operator++();
        
        /**
         * Constructor for a new iterator.
         */
        EdgeIterator(ConformingTetrahedralMesh& rMesh, unsigned elemIndex);
        
    private:
        /** Keep track of what edges have been visited */
        std::set<std::set<unsigned> > mEdgesVisited;
    
        ConformingTetrahedralMesh& mrMesh;
        
        unsigned mElemIndex;
        unsigned mNodeALocalIndex;
        unsigned mNodeBLocalIndex;
        unsigned mCellIndex;
        unsigned mNodeIndex;
    };

    /**
     * @return iterator pointing to the first edge (ie connection between 2 nodes) of the mesh
     */
    EdgeIterator EdgesBegin();
    
    /**
     * @return iterator pointing to one past the last edge (ie connection between 2 nodes) 
     * of the mesh
     */
    EdgeIterator EdgesEnd();
};


/////////////////////////////////////////////////////////////////////////////////////
//   IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConformingTetrahedralMesh()
{
    Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConformingTetrahedralMesh(unsigned numElements)
{
    Clear();
    mElements.reserve(numElements);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConformingTetrahedralMesh(std::vector<Node<SPACE_DIM> *> nodes)
  //: mNodes(nodes)
{
    Clear();
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<SPACE_DIM>* temp_node = nodes[index];
        mNodes.push_back(temp_node);
    }
    mAddedNodes = true;
    NodeMap node_map(nodes.size());
    ReMesh(node_map);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    bool cullInternalFaces)
{
    if(ELEMENT_DIM==1)
    {
        cullInternalFaces = true;
    }

    // Record number of corner nodes
    unsigned num_nodes = rMeshReader.GetNumNodes();
    
    // Reserve memory for nodes, so we don't have problems with pointers stored in
    // elements becoming invalid.
    mNodes.reserve(num_nodes);
    
    rMeshReader.Reset();
    
    //typename std::map<std::pair<unsigned,unsigned>,unsigned>::const_iterator iterator;
    //std::map<std::pair<unsigned,unsigned>,unsigned> internal_nodes_map;
    
    // Add corner nodes
    std::vector<double> coords;
    for (unsigned i=0; i < num_nodes; i++)
    {
        coords = rMeshReader.GetNextNode();
        mNodes.push_back(new Node<SPACE_DIM>(i, coords, false));
    }
    
    //unsigned new_node_index = mNumCornerNodes;
    
    rMeshReader.Reset();
    // Add elements
    //new_node_index = mNumCornerNodes;
    mElements.reserve(rMeshReader.GetNumElements());
    
    for (unsigned element_index=0; element_index < (unsigned) rMeshReader.GetNumElements(); element_index++)
    {
        std::vector<unsigned> node_indices = rMeshReader.GetNextElement();
        std::vector<Node<SPACE_DIM>*> nodes;
        unsigned nodes_size = node_indices.size();
        for (unsigned j=0; j<nodes_size; j++)
        {
            assert(node_indices[j] <  mNodes.size());
            nodes.push_back(mNodes[node_indices[j]]);
        }
        
        mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes));
    }

    
    // Add boundary elements & nodes
    unsigned actual_face_index=0;
    for (unsigned face_index=0; face_index<(unsigned)rMeshReader.GetNumFaces(); face_index++)
    {
        std::vector<unsigned> node_indices = rMeshReader.GetNextFace();
        bool is_boundary_face = true;

        
        // Determine if this is a boundary face
        std::set<unsigned> containing_element_indices; // Elements that contain this face
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned node_index=0; node_index<node_indices.size(); node_index++)
        {
            assert(node_indices[node_index] <  mNodes.size());
            // Add Node pointer to list for creating an element
            nodes.push_back(mNodes[node_indices[node_index]]);
            
            if(cullInternalFaces)
            {                
                // Work out what elements contain this face, by taking the intersection
                // of the sets of elements containing each node in the face.
                if (node_index == 0)
                {
                    containing_element_indices = nodes[node_index]->rGetContainingElementIndices();
                }
                else
                {
                    std::set<unsigned> temp;
                    std::set_intersection(nodes[node_index]->rGetContainingElementIndices().begin(),
                                          nodes[node_index]->rGetContainingElementIndices().end(),
                                          containing_element_indices.begin(), containing_element_indices.end(),
                                          std::inserter(temp, temp.begin()));
                    containing_element_indices = temp;
                }
            }
        }
            
        if(cullInternalFaces)
        {
            //If the following assertion is thrown, it means that the .edge/.face file does not
            //match the .ele file -- they were generated at separate times.  Simply remove the internal
            //edges/faces by hand.
            assert(containing_element_indices.size() != 0);
        
            if(containing_element_indices.size() > 1)
            {
                is_boundary_face = false;
            }
        }

        if (is_boundary_face)
        {
            // This is a boundary face
            // Ensure all its nodes are marked as boundary nodes
            for (unsigned j=0; j<nodes.size(); j++)
            {
                if (!nodes[j]->IsBoundaryNode())
                {
                    nodes[j]->SetAsBoundaryNode();
                    mBoundaryNodes.push_back(nodes[j]);
                }
                //Register the index that this bounday element will have
                //with the node
                nodes[j]->AddBoundaryElement(actual_face_index);
            }
            
            // The added elements will be deleted in our destructor
            mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(actual_face_index, nodes));
            actual_face_index++;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~ConformingTetrahedralMesh()
{
    // Iterate over nodes and free the memory
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        delete mNodes[i];
    }
    // Iterate over elements and free the memory
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    // Iterate over boundary elements and free the memory
    for (unsigned i=0; i<mBoundaryElements.size(); i++)
    {
        delete mBoundaryElements[i];
    }
    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM> *pNewNode)
{

    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(mNodes.size());
        mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index=mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete mNodes[index];
        mNodes[index] = pNewNode;
    }
    mAddedNodes = true;
    return pNewNode->GetIndex();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile)
{
    mNodesPerProcessor.clear();

    std::ifstream file_stream(nodesPerProcessorFile.c_str());
    if(file_stream.is_open())
    {  
        while(file_stream) 
        {
            unsigned nodes_per_processor;
            file_stream >> nodes_per_processor;

            if(file_stream)
            {
                mNodesPerProcessor.push_back(nodes_per_processor);
            }
        }
    }
    else
    {
        EXCEPTION("Unable to read nodes per processor file "+nodesPerProcessorFile);
    }

    unsigned sum = 0;
    for(unsigned i=0; i<mNodesPerProcessor.size(); i++)
    {
        sum += mNodesPerProcessor[i];
    }
    
    if(sum != GetNumNodes())
    {
        std::stringstream string_stream;
        string_stream << "Sum of nodes per processor, " << sum 
                     << ", not equal to number of nodes in mesh, " << GetNumNodes();
        EXCEPTION(string_stream.str());
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned>& ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodesPerProcessor()
{
    return mNodesPerProcessor;
}

/**
 * Get a node reference from the mesh.
 *
 * Note that this may become invalid if nodes are subsequently added to the mesh.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index)
{
    assert(index < mNodes.size());
    return (mNodes[index]);
}

/// Returns the number of nodes that are actually in use
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNodes.size() - mDeletedNodeIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes() const
{
    return mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mElements.size() - mDeletedElementIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements()
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryNodes()
{
    return mBoundaryNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements()
{
    return mBoundaryElements.size() - mDeletedBoundaryElementIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllBoundaryElements()
{
    return mBoundaryElements.size();
}






template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RescaleMeshFromBoundaryNode(ChastePoint<1> updatedPoint, unsigned boundaryNodeIndex)
{
    assert(GetNode(boundaryNodeIndex)->IsBoundaryNode());
    double scaleFactor = updatedPoint[0] / GetNode(boundaryNodeIndex)->GetPoint()[0];
    double temp;
    for (unsigned i=0; i < boundaryNodeIndex+1; i++)
    {
        temp = scaleFactor * mNodes[i]->GetPoint()[0];
        ChastePoint<1> newPoint(temp);
        mNodes[i]->SetPoint(newPoint);
    }
    RefreshMesh();
}

/** 
 *  SetNode moves the node with a particular index to a new point in space and
  * verifies that the signed areas of the supporting Elements are positive
  * @param index is the index of the node to be moved
  * @param point is the new target location of the node
  * @param concreteMove is set to false if we want to skip the signed area tests
  */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned index,
        ChastePoint<SPACE_DIM> point,
        bool concreteMove)
{
    mNodes[index]->SetPoint(point);
    if (concreteMove)
    {
        for (typename Node<SPACE_DIM>::ContainingElementIterator it = mNodes[index]->ContainingElementsBegin();
             it != mNodes[index]->ContainingElementsEnd();
             ++it)
        {
            try
            {
                GetElement(*it)->RefreshJacobianDeterminant();
            }
            catch (Exception e)
            {
                if (ELEMENT_DIM == SPACE_DIM)
                {
                    EXCEPTION("Moving node caused an element to have a non-positive Jacobian determinant");
                }
                else
                {
                    EXCEPTION("Moving node caused an subspace element to change direction");
                }
            }
        }
        for (typename Node<SPACE_DIM>::ContainingBoundaryElementIterator it = mNodes[index]->ContainingBoundaryElementsBegin();
             it != mNodes[index]->ContainingBoundaryElementsEnd();
             ++it)
        {
            try
            {
                GetBoundaryElement(*it)->RefreshJacobianDeterminant();
            }
            catch (Exception e)
            {
                EXCEPTION("Moving node caused a boundary element to have a non-positive Jacobian determinant");
            }
        }
    }
}

/**
 * DeleteNode deletes a node from the mesh by finding an appropriate neighbour node
 * to merge it with.
 *
 * @param index is the index of the node to be deleted
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DeleteNode(unsigned index)
{
    if (mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to delete a deleted node");
    }
    unsigned target_index;
    bool found_target=false;
    for (typename Node<SPACE_DIM>::ContainingElementIterator it = mNodes[index]->ContainingElementsBegin();
         !found_target && it != mNodes[index]->ContainingElementsEnd();
         ++it)
    {
        Element <ELEMENT_DIM,SPACE_DIM> *p_element = GetElement(*it);
        for (unsigned i=0; i<=ELEMENT_DIM && !found_target; i++)
        {
            target_index = p_element->GetNodeGlobalIndex(i);
            try
            {
                MoveMergeNode(index, target_index, false);
                found_target = true;
            }
            catch (Exception e)
            {
                // Just try the next node
            }
        }
    }
    if (!found_target)
    {
        EXCEPTION("Failure to delete node");
    }
    
    MoveMergeNode(index, target_index);
}

/**
 * This marks a node as deleted. Note that it DOES NOT deal with the 
 * associated elements and therefore should only be called immediately prior
 * to a ReMesh() being called. (Thus saves work compared to DeleteNode() 
 * function and does not MoveMerge the node and elements).
 * 
 * @param index The index of the node to delete
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DeleteNodePriorToReMesh(unsigned index)
{
#define COVERAGE_IGNORE
    // A ReMesh can only happen in 2D or 3D so
    assert(SPACE_DIM==2 || SPACE_DIM==3);
#undef COVERAGE_IGNORE
    mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);
}

/**
 * MoveMergeNode moves one node to another (i.e. merges the nodes), refreshing/deleting elements as
 * appropriate.
 *
 * @param index is the index of the node to be moved
 * @param targetIndex is the index of the node to move to
 * @param concreteMove can be set to false if you just want to check whether this will work.
 *     Set it to true if you're doing the merger for real, in order to do all the bookkeeping.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MoveMergeNode(unsigned index,
        unsigned targetIndex,
        bool concreteMove)
{

    if (mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to move a deleted node");
    }
    
    if (index == targetIndex)
    {
        EXCEPTION("Trying to merge a node with itself");
    }
    if (mNodes[index]->IsBoundaryNode())
    {
        if (!mNodes[targetIndex]->IsBoundaryNode())
        {
            EXCEPTION("A boundary node can only be moved on to another boundary node");
        }
    }
    std::set<unsigned> unshared_element_indices;
    std::set_difference(mNodes[index]->rGetContainingElementIndices().begin(),
                        mNodes[index]->rGetContainingElementIndices().end(),
                        mNodes[targetIndex]->rGetContainingElementIndices().begin(),
                        mNodes[targetIndex]->rGetContainingElementIndices().end(),
                        std::inserter(unshared_element_indices, unshared_element_indices.begin()));
                        
                        
    if (unshared_element_indices.size() == mNodes[index]->rGetContainingElementIndices().size())
    {
        EXCEPTION("These nodes cannot be merged since they are not neighbours");
    }
    
    std::set<unsigned> unshared_boundary_element_indices;
    std::set_difference(mNodes[index]->rGetContainingBoundaryElementIndices().begin(),
                        mNodes[index]->rGetContainingBoundaryElementIndices().end(),
                        mNodes[targetIndex]->rGetContainingBoundaryElementIndices().begin(),
                        mNodes[targetIndex]->rGetContainingBoundaryElementIndices().end(),
                        std::inserter(unshared_boundary_element_indices, unshared_boundary_element_indices.begin()));
                        
                        
    if (mNodes[index]->IsBoundaryNode())
    {
        if (unshared_boundary_element_indices.size()
            == mNodes[index]->rGetContainingBoundaryElementIndices().size())
        {
            //May be redundant (only thrown in 1D tests)
            EXCEPTION("These nodes cannot be merged since they are not neighbours on the boundary");
        }
    }
    
    mNodes[index]->rGetModifiableLocation() = mNodes[targetIndex]->rGetLocation();
    
    for (std::set<unsigned>::const_iterator element_iter=unshared_element_indices.begin();
             element_iter != unshared_element_indices.end();
             element_iter++)
    {
        try
        {
        
            GetElement(*element_iter)->RefreshJacobianDeterminant(concreteMove);
            if (concreteMove)
            {
                GetElement(*element_iter)->ReplaceNode(mNodes[index], mNodes[targetIndex]);
            }
            
        }
        catch (Exception e)
        {
            EXCEPTION("Moving node caused an element to have a non-positive Jacobian determinant");
        }
    }
    
    for (std::set<unsigned>::const_iterator boundary_element_iter=
                 unshared_boundary_element_indices.begin();
             boundary_element_iter != unshared_boundary_element_indices.end();
             boundary_element_iter++)
    {
    
        GetBoundaryElement(*boundary_element_iter)->RefreshJacobianDeterminant(concreteMove);
        if (concreteMove)
        {
            GetBoundaryElement(*boundary_element_iter)->ReplaceNode(mNodes[index], mNodes[targetIndex]);
        }
    }
        
    std::set<unsigned> shared_element_indices;
    std::set_intersection(mNodes[index]->rGetContainingElementIndices().begin(),
                          mNodes[index]->rGetContainingElementIndices().end(),
                          mNodes[targetIndex]->rGetContainingElementIndices().begin(),
                          mNodes[targetIndex]->rGetContainingElementIndices().end(),
                          std::inserter(shared_element_indices, shared_element_indices.begin()));
    
    for (std::set<unsigned>::const_iterator element_iter=shared_element_indices.begin();
             element_iter != shared_element_indices.end();
             element_iter++)
    {
        if (concreteMove)
        {
            GetElement(*element_iter)->MarkAsDeleted();
            mDeletedElementIndices.push_back(*element_iter);
        }
        else
        {
            GetElement(*element_iter)->ZeroJacobianDeterminant();
        }
    }
        
        
    std::set<unsigned> shared_boundary_element_indices;
    std::set_intersection(mNodes[index]->rGetContainingBoundaryElementIndices().begin(),
                          mNodes[index]->rGetContainingBoundaryElementIndices().end(),
                          mNodes[targetIndex]->rGetContainingBoundaryElementIndices().begin(),
                          mNodes[targetIndex]->rGetContainingBoundaryElementIndices().end(),
                          std::inserter(shared_boundary_element_indices, shared_boundary_element_indices.begin()));
   
    for (std::set<unsigned>::const_iterator boundary_element_iter=shared_boundary_element_indices.begin();
             boundary_element_iter != shared_boundary_element_indices.end();
             boundary_element_iter++)
    {
        if (concreteMove)
        {
            GetBoundaryElement(*boundary_element_iter)->MarkAsDeleted();
            mDeletedBoundaryElementIndices.push_back(*boundary_element_iter);
        }
        else
        {
            GetBoundaryElement(*boundary_element_iter)->ZeroJacobianDeterminant();
            GetBoundaryElement(*boundary_element_iter)->ZeroWeightedDirection();
        }
    }
        
    if (concreteMove)
    {
        mNodes[index]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(index);
    }
}


/**
 * This method allows the mesh properties to be re-calculated after one
 * or more node have been moved.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
    for (unsigned i=0; i<mElements.size();i++)
    {
        if (!mElements[i]->IsDeleted())
        {
            mElements[i]->RefreshJacobianDeterminant();
        }
    }
    
    //Refresh each boundary element
    for (unsigned i=0; i<mBoundaryElements.size();i++)
    {
        if (!mBoundaryElements[i]->IsDeleted())
        {
            try
            {
                mBoundaryElements[i]->RefreshJacobianDeterminant();
            }
            catch (Exception e)
            {
                //Since we may have rotated the mesh, it's okay for normals to swing round
            }
        }
    }
    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefineElement(
    Element<ELEMENT_DIM,SPACE_DIM>* pElement,
    ChastePoint<SPACE_DIM> point)
{

    //Check that the point is in the element
    if (pElement->IncludesPoint(point, true) == false)
    {
        EXCEPTION("RefineElement could not be started (point is not in element)");
    }
    
    // Add a new node from the point that is passed to RefineElement
    unsigned new_node_index = AddNode(new Node<SPACE_DIM>(0, point.rGetLocation()));
    // Note: the first argument is the index of the node, which is going to be
    //       overriden by AddNode, so it can safely be ignored
    
    //This loop constructs the extra elements which are going to fill the space
    for (unsigned i = 0; i < ELEMENT_DIM; i++)
    {
    
        // First, make a copy of the current element making sure we update its index
        unsigned new_elt_index;
        if (mDeletedElementIndices.empty())
        {
            new_elt_index = mElements.size();
        }
        else
        {
            new_elt_index = mDeletedElementIndices.back();
            mDeletedElementIndices.pop_back();
        }
        
        Element<ELEMENT_DIM,SPACE_DIM>* p_new_element=
            new Element<ELEMENT_DIM,SPACE_DIM>(*pElement, new_elt_index);
            
        // Second, update the node in the element with the new one
        p_new_element->UpdateNode(ELEMENT_DIM-1-i, mNodes[new_node_index]);
        
        p_new_element->RefreshJacobianDeterminant();
        
        // Third, add the new element to the set
        if ((unsigned) new_elt_index == mElements.size())
        {
            mElements.push_back(p_new_element);
        }
        else
        {
            delete mElements[new_elt_index];
            mElements[new_elt_index] = p_new_element;
        }
        
    }
    
    // Lastly, update the last node in the element to be refined
    pElement->UpdateNode(ELEMENT_DIM, mNodes[new_node_index]);
    pElement->RefreshJacobianDeterminant();
    
    return new_node_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateVolume()
{
    double mesh_volume = 0.0;
    
    ElementIterator it = GetElementIteratorBegin();
    
    while (it != GetElementIteratorEnd())
    {
        mesh_volume += (*it)->GetVolume();
        it++;
    }
    
    return mesh_volume;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateSurfaceArea()
{
    //ELEMENT_DIM-1 is the dimension of the boundary element
    assert (ELEMENT_DIM>=1);
    const unsigned bound_element_dim=ELEMENT_DIM-1;
    assert(bound_element_dim < 3);
    if ( bound_element_dim == 0)
    {
        return 0.0;
    }
    
    double mesh_surface = 0.0;
    BoundaryElementIterator it = GetBoundaryElementIteratorBegin();
    
    while (it != GetBoundaryElementIteratorEnd())
    {
        mesh_surface += (*it)->GetJacobianDeterminant();
        it++;
    }
    
    if ( bound_element_dim == 2)
    {
        mesh_surface /= 2.0;
    }
    
    return mesh_surface;
}



/**
 * Scale the mesh.
 * @param xFactor is the scale in the x-direction,
 * @param yFactor is the scale in the y-direction,
 * @param zFactor is the scale in the z-direction
 **/
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Scale(
    const double xScale,
    const double yScale,
    const double zScale)
{
    unsigned num_nodes=GetNumAllNodes();
    
    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = mNodes[i]->rGetModifiableLocation();
        if (SPACE_DIM>=3)
        {
            r_location[2] *= zScale;
        }
        if (SPACE_DIM>=2)
        {
            r_location[1] *= yScale;
        }
        r_location[0] *= xScale;
    }
    
    RefreshMesh();
}

/**
 * Translate the mesh.
 * @param xMovement is the x-displacement,
 * @param yMovement is the y-displacement,
 * @param zMovement is the z-displacement,
 **/
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Translate(
    const double xMovement,
    const double yMovement,
    const double zMovement)
{
    c_vector<double , SPACE_DIM> displacement;
    
    switch (SPACE_DIM)
    {
        case 3:
            displacement[2]=zMovement;
        case 2:
            displacement[1]=yMovement;
        case 1:
            displacement[0]=xMovement;
    }
    
    Translate(displacement);
}


/**
 * Translate mesh using the BOOST ublas library - this is the method that actually does the work
 * @param transVec is a translation vector of the correct size
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Translate(c_vector<double, SPACE_DIM> transVec)
{
    unsigned num_nodes=GetNumAllNodes();
    
    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = mNodes[i]->rGetModifiableLocation();
        r_location += transVec;
    }
    
    RefreshMesh();
}




/**
 * Do a general mesh rotation with an +ve determinant orthonormal rotation_matrix
 * This is the method that actually does the work
 * @param rotation_matrix is a Ublas rotation matrix of the correct form
 **/
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(
    c_matrix<double , SPACE_DIM, SPACE_DIM> rotation_matrix)
{
    unsigned num_nodes=GetNumAllNodes();
    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = mNodes[i]->rGetModifiableLocation();
        r_location = prod(rotation_matrix, r_location);
    }
    
    RefreshMesh();
}

/**
* Do an angle axis rotation
* @param axis is the axis of rotation (does not need to be normalised)
* @param angle is the angle in radians
 **/
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(c_vector<double,3> axis, double angle)
{
    assert(SPACE_DIM == 3);
    double norm = norm_2(axis);
    c_vector<double,3> unit_axis=axis/norm;
    
    c_matrix<double, SPACE_DIM,SPACE_DIM> rotation_matrix;
    
    double c = cos(angle);
    double s = sin(angle);
    
    rotation_matrix(0,0) = unit_axis(0)*unit_axis(0)+c*(1-unit_axis(0)*unit_axis(0));
    rotation_matrix(0,1) = unit_axis(0)*unit_axis(1)*(1-c) - unit_axis(2)*s;
    rotation_matrix(1,0) = unit_axis(0)*unit_axis(1)*(1-c) + unit_axis(2)*s;
    rotation_matrix(1,1) = unit_axis(1)*unit_axis(1)+c*(1-unit_axis(1)*unit_axis(1));
    rotation_matrix(0,2) = unit_axis(0)*unit_axis(2)*(1-c)+unit_axis(1)*s;
    rotation_matrix(1,2) = unit_axis(1)*unit_axis(2)*(1-c)-unit_axis(0)*s;
    rotation_matrix(2,0) = unit_axis(0)*unit_axis(2)*(1-c)-unit_axis(1)*s;
    rotation_matrix(2,1) = unit_axis(1)*unit_axis(2)*(1-c)+unit_axis(0)*s;
    rotation_matrix(2,2) = unit_axis(2)*unit_axis(2)+c*(1-unit_axis(2)*unit_axis(2));
    
    Rotate(rotation_matrix);
}


/**
 * Rotate the mesh about the x-axis
 * @param theta is the angle in radians
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateX(const double theta)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("This rotation is only valid in 3D");
    }
    c_matrix<double , SPACE_DIM, SPACE_DIM> x_rotation_matrix=identity_matrix<double>(SPACE_DIM);
    
    x_rotation_matrix(1,1) = cos(theta);
    x_rotation_matrix(1,2) = sin(theta);
    x_rotation_matrix(2,1) = -sin(theta);
    x_rotation_matrix(2,2) = cos(theta);
    Rotate(x_rotation_matrix);
}


/**
 * Rotate the mesh about the y-axis
 * @param theta is the angle in radians
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateY(const double theta)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("This rotation is only valid in 3D");
    }
    c_matrix<double , SPACE_DIM, SPACE_DIM> y_rotation_matrix=identity_matrix<double>(SPACE_DIM);
    
    y_rotation_matrix(0,0) = cos(theta);
    y_rotation_matrix(0,2) = -sin(theta);
    y_rotation_matrix(2,0) = sin(theta);
    y_rotation_matrix(2,2) = cos(theta);
    
    
    Rotate(y_rotation_matrix);
}

/**
 * Rotate the mesh about the z-axis
 * @param theta is the angle in radians
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateZ(const double theta)
{
    if (SPACE_DIM < 2)
    {
        EXCEPTION("This rotation is not valid in less than 2D");
    }
    c_matrix<double , SPACE_DIM, SPACE_DIM> z_rotation_matrix=identity_matrix<double>(SPACE_DIM);
    
    
    z_rotation_matrix(0,0) = cos(theta);
    z_rotation_matrix(0,1) = sin(theta);
    z_rotation_matrix(1,0) = -sin(theta);
    z_rotation_matrix(1,1) = cos(theta);
    
    Rotate(z_rotation_matrix);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DeleteBoundaryNodeAt(unsigned index)
{
    if (!mNodes[index]->IsBoundaryNode() )
    {
        EXCEPTION(" You may only delete a boundary node ");
    }
    
    mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);
    // Update the boundary node vector
    typename std::vector<Node<SPACE_DIM>*>::iterator b_node_iter
    = std::find(mBoundaryNodes.begin(), mBoundaryNodes.end(), mNodes[index]);
    mBoundaryNodes.erase(b_node_iter);
    
    // Remove boundary elements containing this node
    std::set<unsigned> boundary_element_indices = mNodes[index]->rGetContainingBoundaryElementIndices();
    std::set<unsigned>::const_iterator boundary_element_indices_iterator = boundary_element_indices.begin();
    while (boundary_element_indices_iterator != boundary_element_indices.end())
    {
        BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element = GetBoundaryElement(*boundary_element_indices_iterator);
        p_boundary_element->MarkAsDeleted();
        mDeletedBoundaryElementIndices.push_back(*boundary_element_indices_iterator);
        boundary_element_indices_iterator++;
    }
    
    // Remove elements containing this node
    std::set<unsigned> element_indices = mNodes[index]->rGetContainingElementIndices();
    std::set<unsigned>::const_iterator element_indices_iterator = element_indices.begin();
    while (element_indices_iterator != element_indices.end())
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(*element_indices_iterator);
        for (unsigned i=0 ; i< p_element->GetNumNodes();i++)
        {
            Node<SPACE_DIM>* p_node = p_element->GetNode(i);
            if (!p_node->IsDeleted())
            {
                p_node->SetAsBoundaryNode();
                // Update the boundary node vector
                mBoundaryNodes.push_back(p_node);
            }
        }
        p_element->MarkAsDeleted();
        mDeletedElementIndices.push_back(p_element->GetIndex());
        element_indices_iterator++;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReIndex(NodeMap& map)
{
    assert(!mAddedNodes);
    map.Resize(GetNumAllNodes());
    
    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> live_elements;
    for (unsigned i=0; i<mElements.size(); i++)
    {
        if (mElements[i]->IsDeleted())
        {
            delete mElements[i];
        }
        else
        {
            live_elements.push_back(mElements[i]);
        }
    }
    
    assert (mDeletedElementIndices.size() == mElements.size()-live_elements.size());    
    mDeletedElementIndices.clear();
    mElements = live_elements;
    
    std::vector<Node<SPACE_DIM> *> live_nodes;
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        if (mNodes[i]->IsDeleted())
        {
            delete mNodes[i];
            map.SetDeleted(i);
        }
        else
        {
            live_nodes.push_back(mNodes[i]);
            // the nodes will have their index set to be the index into the live_nodes
            // vector further down
            map.SetNewIndex(i, (unsigned)(live_nodes.size()-1));
        }
    }
    
    assert (mDeletedNodeIndices.size() == mNodes.size()-live_nodes.size());
    mNodes = live_nodes;
    mDeletedNodeIndices.clear();
    
    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> live_boundary_elements;
    for (unsigned i=0; i<mBoundaryElements.size(); i++)
    {
        if (mBoundaryElements[i]->IsDeleted())
        {
            delete mBoundaryElements[i];
        }
        else
        {
            live_boundary_elements.push_back(mBoundaryElements[i]);
        }
    }
    
    assert (mDeletedBoundaryElementIndices.size() == mBoundaryElements.size()-live_boundary_elements.size());
    mBoundaryElements = live_boundary_elements;
    mDeletedBoundaryElementIndices.clear();
 
    for (unsigned i=0; i<mNodes.size();i++)
    {
        mNodes[i]->SetIndex(i);
    }
    for (unsigned i=0; i<mElements.size();i++)
    {
        mElements[i]->ResetIndex(i);
    }
    for (unsigned i=0; i<mBoundaryElements.size();i++)
    {
        mBoundaryElements[i]->ResetIndex(i);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReMeshWithTriangleLibrary(NodeMap& map)
{
    struct triangulateio triangle_input;
    triangle_input.pointlist = (double *) malloc(GetNumNodes() * 2 * sizeof(double));
    triangle_input.numberofpoints = GetNumNodes();
    triangle_input.numberofpointattributes = 0;
    triangle_input.pointmarkerlist = NULL;
    triangle_input.numberofsegments = 0;
    triangle_input.numberofholes = 0;
    triangle_input.numberofregions = 0;
    
    unsigned new_index = 0;
    map.Resize(GetNumAllNodes());        
    for (unsigned i=0; i<GetNumAllNodes(); i++)
    {
        if (mNodes[i]->IsDeleted())
        {
            map.SetDeleted(i);
        }
        else
        {
            map.SetNewIndex(i,new_index);
            triangle_input.pointlist[2*new_index]=mNodes[i]->rGetLocation()[0];
            triangle_input.pointlist[2*new_index + 1]=mNodes[i]->rGetLocation()[1];
            new_index++;
           
        }
    }
    
    //Make structure for output
    struct triangulateio triangle_output;
    triangle_output.pointlist =  NULL;            
    triangle_output.pointattributelist = (double *) NULL;
    triangle_output.pointmarkerlist = (int *) NULL; 
    triangle_output.trianglelist = (int *) NULL;          
    triangle_output.triangleattributelist = (double *) NULL;
    triangle_output.edgelist = (int *) NULL;             
    triangle_output.edgemarkerlist = (int *) NULL;  

    //Library call 
    triangulate((char*)"Qze", &triangle_input, &triangle_output, NULL);
    
    assert(triangle_output.numberofcorners == 3);
    
    //Remove current data
    Clear();
    
    //Construct the nodes
    for (unsigned node_index=0; node_index<(unsigned)triangle_output.numberofpoints; node_index++)
    {
        if (triangle_output.pointmarkerlist[node_index] == 1)
        {
            //Boundary node
            Node<SPACE_DIM> *p_node=new Node<SPACE_DIM>(node_index, true, 
              triangle_output.pointlist[node_index * 2], 
              triangle_output.pointlist[node_index * 2+1]);
            mNodes.push_back(p_node);
            mBoundaryNodes.push_back(p_node);
        }
        else
        {
            mNodes.push_back(new Node<SPACE_DIM>(node_index, false,
              triangle_output.pointlist[node_index * 2], 
              triangle_output.pointlist[node_index * 2+1]));
        }
        
    }
    
    //Construct the elements
    mElements.reserve(triangle_output.numberoftriangles);
    for (unsigned element_index=0; element_index < (unsigned)triangle_output.numberoftriangles; element_index++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned j = 0; j < 3; j++) 
        {
            unsigned global_node_index=triangle_output.trianglelist[element_index*3 + j];
            assert(global_node_index <  mNodes.size());
            nodes.push_back(mNodes[global_node_index]);
        }
        mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes));
    }
    
    //Construct the edges
    //too big mBoundaryElements.reserve(triangle_output.numberoftriangles);
    unsigned next_boundary_element_index=0;
    for (unsigned boundary_element_index=0; boundary_element_index < (unsigned)triangle_output.numberofedges; boundary_element_index++)
    {
        if (triangle_output.edgemarkerlist[boundary_element_index] == 1)
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            for (unsigned j = 0; j < 2; j++) 
            {
                unsigned global_node_index=triangle_output.edgelist[boundary_element_index*2 + j];
                assert(global_node_index <  mNodes.size());
                nodes.push_back(mNodes[global_node_index]);
            }
            mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(next_boundary_element_index++, nodes));
        }
    }
     
      
    
    
    free(triangle_input.pointlist);
    
    free(triangle_output.pointlist);
    free(triangle_output.pointattributelist);
    free(triangle_output.pointmarkerlist);
    free(triangle_output.trianglelist);
    free(triangle_output.triangleattributelist);
    free(triangle_output.edgelist);
    free(triangle_output.edgemarkerlist);
    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(NodeMap& map)
{
    //Make sure that we are in the correct dimension -- this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 || SPACE_DIM==3 );
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE
    // avoid some triangle/tetgen errors:
    // need at least four nodes for tetgen, and at least three for triangle 
    // assert( GetNumNodes() > SPACE_DIM );

    //Make sure the map is big enough
    map.Resize(GetNumAllNodes());
    

//// #554: commented out for the time being as some cancer tests do not pass with
//// this (probably because of hardcoded values), although it seems to work ok
//// and passes all the mesh tests.
//    if (mDeletedNodeIndices.size()==0 && !mAddedNodes)
//    {
//        //If there are no nodes waiting to be deleted and the current mesh is
//        //Voronoi then we don't need to call triangle/tetgen
//        if (CheckVoronoi())
//        {
//            map.ResetToIdentity();
//            return;
//        }
//    }
    std::stringstream pid;
    pid<<getpid();
    
    OutputFileHandler handler("");
    std::string full_name = handler.GetOutputDirectoryFullPath("")+"temp_"+pid.str()+".";
    
    // Only the master process should do IO and call the mesher
    if (handler.IsMaster())
    {
        std::string node_file_name="temp_"+pid.str()+".node";
        {//Scope for node_file
            out_stream node_file=handler.OpenOutputFile(node_file_name);
        
            (*node_file)<<GetNumNodes()<<"\t" << SPACE_DIM << "\t0\t0\n";
            
            unsigned new_index = 0;
            
            for (unsigned i=0; i<GetNumAllNodes(); i++)
            {
                if (mNodes[i]->IsDeleted())
                {
                    map.SetDeleted(i);
                }
                else
                {
                    map.SetNewIndex(i,new_index);
                    new_index++;
                    const c_vector<double, SPACE_DIM> node_loc = mNodes[i]->rGetLocation();
                    (*node_file)<<i<<"\t"<<node_loc[0]<<"\t"<<node_loc[1];
                    if (SPACE_DIM ==3)
                    {
                        (*node_file)<<"\t"<<node_loc[2];
                    }
                    (*node_file)<<"\n";
                }
            }
            
            node_file->close();
            
        }//Scope for node_file
        
        
        std::string binary_name;
        if (SPACE_DIM==2)
        {
            if (sizeof(long)==4)
            {
                //32-bit machine, so use default binary.
                binary_name="triangle";
            }
            else
            {
                binary_name="triangle_64";
            }
        }
        else
        {
            binary_name="tetgen";
        }
        std::string command =   "./bin/"+ binary_name +" -Qe "
                              + full_name + "node";
        
        if (SPACE_DIM == 3)
        {
            //Tetgen's quiet mode isn't as quiet as Triangle's
            command += " > /dev/null";
        }
        int return_value = system(command.c_str());
        
        
        if (return_value != 0)
        {
            EXCEPTION("The triangle/tetgen mesher did not succeed in remeshing.");
        }
    }
    // Wait for the new mesh to be available and communicate its name
#ifndef SPECIAL_SERIAL
    if (!PetscTools::IsSequential())
    {
        char full_name_comm[200];
        strcpy(full_name_comm, full_name.c_str());
        MPI_Bcast(full_name_comm, 200, MPI_CHAR, 0, MPI_COMM_WORLD);
        full_name=full_name_comm;
    }
#endif //SPECIAL_SERIAL
    
    // clear all current data
    Clear();
    
    //Read the new mesh back from file
    TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(full_name+"1");    
    ConstructFromMeshReader(mesh_reader);
    
    // Make sure the file is not deleted before all the processors have read it
#ifndef SPECIAL_SERIAL
    if (!PetscTools::IsSequential())
    {
        MPI_Barrier(PETSC_COMM_WORLD);
    }
#endif //SPECIAL_SERIAL
    
    if (handler.IsMaster())
    {
        std::string remove_command = "rm "+ full_name+"*";
        system(remove_command.c_str());
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh()
{
    NodeMap map(GetNumNodes());
    ReMesh(map);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    RandomNumberGenerator *p_rng=RandomNumberGenerator::Instance();
    
    //Working from the back, each node is swapped with a random node that precedes it in the array
    for (unsigned index=mNodes.size()-1; index>0; index--)
    {
        unsigned  other=p_rng->randMod(index+1); //includes the possibility of rolling "index"
        //Swap index and other
        Node<SPACE_DIM> *temp=mNodes[index];
        mNodes[index]=mNodes[other];
        mNodes[other]=temp;
    }
    
    //Update indices
    for (unsigned index=0; index<mNodes.size(); index++)
    {
        mNodes[index]->SetIndex(index);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes(std::vector<unsigned>& perm)
{
    //Let's not do this if there are any deleted nodes
    assert( GetNumAllNodes() == GetNumNodes());
    
    assert(perm.size() == mNodes.size());   
    
    //Copy the node pointers
    std::vector <Node <SPACE_DIM> *> copy_m_nodes;
    copy_m_nodes.assign(mNodes.begin(), mNodes.end());
        
    for (unsigned i=0;i<mNodes.size();i++)
    {
        assert(perm[i] < mNodes.size());
        mNodes[ perm[i] ] = copy_m_nodes[i];
    }
    
    //Update indices
    for (unsigned index=0; index<mNodes.size(); index++)
    {
        mNodes[index]->SetIndex(index);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodesWithMetisBinaries(unsigned numProcs)
{
    assert( ELEMENT_DIM==2 || ELEMENT_DIM==3 );
    assert( GetNumAllElements() == GetNumElements());
    assert( GetNumAllNodes() == GetNumNodes());
    
    // Open a file for the elements
    OutputFileHandler handler("");
    
    // Filenames
    std::string basename = "metis.mesh";
    std::stringstream output_file;
    output_file << basename << ".npart." << numProcs;
    std::string nodes_per_proc_file = basename + ".nodesperproc";
    
    // Only the master process should do IO and call METIS
    if (handler.IsMaster())
    {
        out_stream metis_file=handler.OpenOutputFile(basename);
        
        (*metis_file)<<GetNumElements()<<"\t";
        if (ELEMENT_DIM==2)
        {
            (*metis_file)<<1<<"\n"; //1 is Metis speak for triangles
        }
        else
        {
            (*metis_file)<<2<<"\n"; //2 is Metis speak for tetrahedra
        }
        
        for (unsigned i=0; i<(unsigned)GetNumElements(); i++)
        {
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                //Note the +1 since Metis wants meshes indexed from 1
                (*metis_file)<<mElements[i]->GetNode(j)->GetIndex() + 1<<"\t";
            }
            (*metis_file)<<"\n";
        }
        metis_file->close();
                
        /*
         *  Call METIS binary to perform the partitioning. 
         *  It will output a file called metis.mesh.npart.numProcs
         */
        std::stringstream permute_command;
        permute_command <<  "./bin/partdmesh "
                        <<  handler.GetOutputDirectoryFullPath("")
                        <<  basename << " "
                        <<  numProcs
                        <<  " > /dev/null";
                        
        system(permute_command.str().c_str());

        /*
         *  Create a file with the number of nodes per partition
         */
        // Make sure it doesn't exist, since values will be appended with >> 
        std::stringstream clear_command;
        clear_command << "rm -f "
                      << handler.GetOutputDirectoryFullPath("")
                      << nodes_per_proc_file
                      << " > /dev/null";
        system(clear_command.str().c_str());

        // Loop over the partition number (i.e. processor number) and count how many nodes                      
        for (unsigned proc_index=0; proc_index<numProcs; proc_index++)
        {
            std::stringstream count_command;
            count_command << "grep " 
                          << proc_index << " "
                          << handler.GetOutputDirectoryFullPath("")
                          << output_file.str()
                          << " | wc -l >> "
                          << handler.GetOutputDirectoryFullPath("")
                          << nodes_per_proc_file; 
                          
            system(count_command.str().c_str());               
        }

    }
    
    // Wait for the permutation to be available
#ifndef SPECIAL_SERIAL
    PetscTools::Barrier();
#endif    

    /*
     *  Read partition file back into a vector.
     */    
    std::vector<unsigned> partition(GetNumNodes());
    std::vector<unsigned> offset(numProcs,0u);
        
    std::ifstream partition_stream;
    std::string full_path = handler.GetOutputDirectoryFullPath("")
                            + output_file.str(); 
    
    partition_stream.open(full_path.c_str());
    assert(partition_stream.is_open()); 

    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        unsigned part_read;
        
        partition_stream >> part_read;
        
        partition[node_index] = part_read;
        for (unsigned proc=part_read+1; proc<numProcs; proc++)
        {            
            offset[proc]++;
        }
    }
    partition_stream.close();
    
    std::cout << offset [ partition[0] ] << std::endl;
    
    /*
     *  Create the permutation vector based on Metis output
     */    
    std::vector<unsigned> permutation(GetNumNodes(), UINT_MAX);
    std::vector<unsigned> count(numProcs,0u);
    
    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        unsigned part = partition[node_index];
        // Permutation defined like: node number "offset[part] + count[part]" will be renumbered as node_index
        permutation [ node_index ] = offset[part] + count[part];
           
        count[part]++;
    }
    
    PermuteNodes(permutation);    
    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CheckVoronoi(Element<ELEMENT_DIM, SPACE_DIM> *pElement, double maxPenetration)
{
    assert (ELEMENT_DIM == SPACE_DIM);
    unsigned num_nodes = pElement->GetNumNodes();
    std::set<unsigned> neighbouring_elements_indices;
    std::set< Element<ELEMENT_DIM,SPACE_DIM> *> neighbouring_elements;
    std::set<unsigned> neighbouring_nodes_indices;
    
    //Form a set of neighbouring elements via the nodes
    for (unsigned i = 0 ; i < num_nodes; i++)
    {
        Node<SPACE_DIM>* node = pElement->GetNode(i);
        neighbouring_elements_indices = node->rGetContainingElementIndices();
        ///\todo Should use a set union operation here
        for (std::set<unsigned>::const_iterator it = neighbouring_elements_indices.begin();
                 it != neighbouring_elements_indices.end(); ++it)
            {
                neighbouring_elements.insert(GetElement(*it));
            }
    }
    neighbouring_elements.erase(pElement);
    
    //For each neighbouring element find the supporting nodes
    typedef typename std::set<Element<ELEMENT_DIM,SPACE_DIM> *>::const_iterator ElementIterator;
    
    for (ElementIterator it = neighbouring_elements.begin();
         it != neighbouring_elements.end(); ++it)
    {
        for (unsigned i = 0 ; i < num_nodes; i++)
        {
            neighbouring_nodes_indices.insert((*it)->GetNodeGlobalIndex(i));
        }
    }
    //Remove the nodes that support this element
    for (unsigned i = 0 ; i < num_nodes; i++)
    {
        neighbouring_nodes_indices.erase(pElement->GetNodeGlobalIndex(i));
    }
    
    //Get the circumsphere information
    c_vector <double, ELEMENT_DIM+1> this_circum_centre;
    this_circum_centre = pElement->CalculateCircumsphere();
    
    //Copy the actualy circumcentre into a smaller vector
    c_vector <double, ELEMENT_DIM> circum_centre;
    for (unsigned i=0;i<ELEMENT_DIM;i++)
    {
        circum_centre[i]=this_circum_centre[i];
    }
    
    for (std::set<unsigned>::const_iterator it = neighbouring_nodes_indices.begin();
             it != neighbouring_nodes_indices.end(); ++it)
    {
        c_vector <double, ELEMENT_DIM> node_location = GetNode(*it)->rGetLocation();
        
        // Calculate vector from circumcenter to node
        node_location -= circum_centre;
        // This is to calculate the squared distance betweeen them
        double squared_distance = inner_prod(node_location, node_location);
        
        // If the squared idstance is less than the elements circum-radius(squared),
        // then the voronoi property is violated.
        
        if (squared_distance < this_circum_centre[ELEMENT_DIM])
        {
            // We know the node is inside the circumsphere, but we don't know how far
            double radius = sqrt(this_circum_centre[ELEMENT_DIM]);
            double distance = radius - sqrt(squared_distance);
            
            // If the node penetration is greater than supplied maximum penetration factor
            if (distance/radius > maxPenetration)
            {
                return false;
            }
        }
    }
    return true;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CheckVoronoi(double maxPenetration)
{
    // Looping through all the elements in the mesh
    for (unsigned i=0; i < mElements.size();i++)
    {
        // Check if the element is not deleted
        if (!mElements[i]->IsDeleted())
        {
            // Checking the Voronoi of the Element
            if (CheckVoronoi(mElements[i], maxPenetration) == false)
            {
                return false;
            }
        }
    }
    return true;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructLinearMesh(unsigned width)
{
    assert(SPACE_DIM == 1);
    assert(ELEMENT_DIM == 1);
    
    for (unsigned node_index=0; node_index<=width; node_index++)
    {
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, node_index==0 || node_index==width, node_index);
        mNodes.push_back(p_node); // create node
        if (node_index==0) // create left boundary node and boundary element
        {
            mBoundaryNodes.push_back(p_node);
            mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(0, p_node) );
        }
        if (node_index==width) // create right boundary node and boundary element
        {
            mBoundaryNodes.push_back(p_node);
            mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(1, p_node) );
        }
        if (node_index>0) // create element
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            nodes.push_back(mNodes[node_index-1]);
            nodes.push_back(mNodes[node_index]);
            mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(node_index-1, nodes) );
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRectangularMesh(unsigned width, unsigned height, bool stagger)
{
    assert(SPACE_DIM == 2);
    assert(ELEMENT_DIM == 2);
    
    //Construct the nodes
    unsigned node_index=0;
    for (int j=(int)height;j>=0;j--) //j must be signed for this loop to terminate
    {
        for (unsigned i=0;i<width+1;i++)
        {
            bool is_boundary=false;
            if (i==0 || j==0 || i==width || j==(int)height)
            {
                is_boundary=true;
            }
            Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index++, is_boundary, i, j);
            mNodes.push_back(p_node);
            if(is_boundary)
            {
                mBoundaryNodes.push_back(p_node);
            }
        }
    }
    
    //Construct the boundary elements
    unsigned belem_index=0;
    //Top
    for (unsigned i=0;i<width;i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(mNodes[i]);
        nodes.push_back(mNodes[i+1]);
        mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Right
    for (unsigned i=1;i<height+1;i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(mNodes[(width+1)*i-1]);
        nodes.push_back(mNodes[(width+1)*(i+1)-1]);
        mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Bottom
    for (unsigned i=0;i<width;i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(mNodes[height*(width+1)+i+1]);
        nodes.push_back(mNodes[height*(width+1)+i]);
        mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Left
    for (unsigned i=0;i<height;i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(mNodes[(width+1)*(i+1)]);
        nodes.push_back(mNodes[(width+1)*(i)]);
        mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    
    //Construct the elements
    unsigned elem_index=0;
    for (unsigned j=0;j<height;j++)
    {
        for (unsigned i=0; i<width; i++)
        {
            unsigned parity=(i+j)%2;
            std::vector<Node<SPACE_DIM>*> upper_nodes;
            upper_nodes.push_back(mNodes[j*(width+1)+i]);
            upper_nodes.push_back(mNodes[j*(width+1)+i+1]);
            if (stagger==false  || parity == 0)
            {
                upper_nodes.push_back(mNodes[(j+1)*(width+1)+i+1]);
            }
            else
            {
                upper_nodes.push_back(mNodes[(j+1)*(width+1)+i]);
            }
            mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,upper_nodes));
            std::vector<Node<SPACE_DIM>*> lower_nodes;
            lower_nodes.push_back(mNodes[(j+1)*(width+1)+i+1]);
            lower_nodes.push_back(mNodes[(j+1)*(width+1)+i]);
            if (stagger==false  ||parity == 0)
            {
                lower_nodes.push_back(mNodes[j*(width+1)+i]);
            }
            else
            {
                lower_nodes.push_back(mNodes[j*(width+1)+i+1]);
            }
            mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,lower_nodes));
        }
    }
    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndex(ChastePoint<SPACE_DIM> testPoint, bool strict, std::set<unsigned> testElements)
{
    for(std::set<unsigned>::iterator iter=testElements.begin(); iter!=testElements.end(); iter++)
    {
        assert(*iter<GetNumElements());
        ///\todo What if the element is deleted?
        if (mElements[*iter]->IncludesPoint(testPoint, strict))
        {
            return *iter;
        }
    }        
    
    ///\todo This ought to return a set of all elements that contain the point (if the point is a node in the mesh then it's contained in multiple elements)
    ///\todo Polling every element is unnecessary.  We ought to start from a likely place and hill climb
    for (unsigned i=0; i<mElements.size(); i++)
    {
        ///\todo What if the element is deleted?
        if (mElements[i]->IncludesPoint(testPoint, strict))
        {
            return i;
        }
    }
    
    //If it's in none of the elements, then throw
    EXCEPTION("Point is not in mesh");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNearestElementIndex(ChastePoint<SPACE_DIM> testPoint)
{
    ///\todo This ought to return a set of all elements that contain the point (if the point is a node in the mesh then it's contained in multiple elements)
    ///\todo Polling every element is unnecessary.  We ought to start from a likely place and hill climb
    
    double max_min_weight=-INFINITY;
    unsigned closest_index=0;
    for (unsigned i=0; i < mElements.size();i++)
    {
        ///\todo What if the element is deleted?
        c_vector<double, ELEMENT_DIM+1> weight=mElements[i]->CalculateInterpolationWeights(testPoint);
        double min_weight=1.0;
        for (unsigned j=0; j<=ELEMENT_DIM; j++)
        {
            if (weight[j]<min_weight)
            {
                min_weight=weight[j];
            }
        }
        if (min_weight > max_min_weight)
        {
            max_min_weight = min_weight;
            closest_index=i;
        }
        
    }
    return closest_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndices(ChastePoint<SPACE_DIM> testPoint)
{
    std::vector<unsigned> element_indices;
    for (unsigned i=0; i < mElements.size();i++)
    {
        ///\todo What if the element is deleted?
        if (mElements[i]->IncludesPoint(testPoint))
        {
            element_indices.push_back(i);
        }
    }
    return element_indices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
{
    assert(hi>=lo);
    for (unsigned element_index=0; element_index<mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element=mElements[element_index];
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructCuboid(unsigned width,
        unsigned height,
        unsigned depth,
        bool stagger)
{
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 3);
    //Construct the nodes
    
    unsigned node_index=0;
    for (unsigned k=0;k<depth+1;k++)
    {
        for (unsigned j=0;j<height+1;j++)
        {
            for (unsigned i=0;i<width+1;i++)
            {
                bool is_boundary=false;
                if (i==0 || j==0 || k==0 || i==width || j==height || k==depth)
                {
                    is_boundary=true;
                }
                
                Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index++, is_boundary, i, j, k);

                mNodes.push_back(p_node);
                if(is_boundary)
                {
                    mBoundaryNodes.push_back(p_node);
                }
            }
        }
    }
    
    // Construct the elements
    
    unsigned elem_index=0;
    unsigned belem_index=0;
    unsigned element_nodes[4][6][4] = {{{0, 1, 5, 7}, {0, 1, 3, 7},
                                        {0, 2, 3, 7}, {0, 2, 6, 7},
                                        {0, 4, 6, 7}, {0, 4, 5, 7}},
                                       {{1, 0, 2, 6}, {1, 0, 4, 6},
                                        {1, 5, 4, 6}, {1, 5, 7, 6},
                                        {1, 3, 2, 6}, {1, 3, 7, 6}},
                                       {{2, 0, 1, 5}, {2, 0, 4, 5},
                                        {2, 3, 1, 5}, {2, 3, 7, 5},
                                        {2, 6, 4, 5}, {2, 6, 7, 5}},
                                       {{3, 1, 0, 4}, {3, 1, 5, 4},
                                        {3, 2, 0, 4}, {3, 2, 6, 4},
                                        {3, 7, 5, 4}, {3, 7, 6, 4}}};     
                                    
    std::vector<Node<SPACE_DIM>*> tetrahedra_nodes;
    
    for (unsigned k=0;k<depth;k++)
    {
        for (unsigned j=0;j<height;j++)
        {
            for (unsigned i=0;i<width;i++)
            {
                // Compute the nodes' index
                unsigned global_node_indices[8];
                unsigned local_node_index = 0;
                
                for (unsigned z = 0; z < 2; z++)
                {
                    for (unsigned y = 0; y < 2; y++)
                    {
                        for (unsigned x = 0; x < 2; x++)
                        {
                            global_node_indices[local_node_index] = i+x+(width+1)*(j+y+(height+1)*(k+z));
                            
                            local_node_index++;
                        }
                    }
                }
                
                for (unsigned m = 0; m < 6; m++)
                {
                    // Tetrahedra #m
                    
                    tetrahedra_nodes.clear();
                                        
                    for (unsigned n = 0; n < 4; n++)
                    {
                        if (stagger)
                        {
                            if (i%2==0)
                            {
                                if (j%2==0)
                                {
                                    if (k%2==0)
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[0][m][n]]]);    
                                    }
                                    else
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[3][m][n]]]);
                                    }
                                }
                                else
                                {
                                    if (k%2==0)
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[2][m][n]]]);
                                    }
                                    else
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[1][m][n]]]);
                                    }
                                }
                            }
                            else
                            {
                               if (j%2==0)
                               {
                                    if (k%2==0)
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[3][m][n]]]);
                                    }
                                    else
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[1][m][n]]]);
                                    }
                               }
                               else
                               {
                                    if (k%2==0)
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[2][m][n]]]);
                                    }
                                    else
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[0][m][n]]]);
                                    } 
                                }
                            }
                        }
                                    
                        else
                        {
                            tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[0][m][n]]]);
                        }
                    }
                    
                    mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++, tetrahedra_nodes));
                }
                
                //Are we at a boundary?
                std::vector<Node<SPACE_DIM>*> triangle_nodes;
                if (i == 0) //low face at x==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[6]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[6]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[4]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (i == width-1) //high face at x=width
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[5]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[3]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == 0) //low face at y==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[5]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[1]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[5]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == height-1) //high face at y=height
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[3]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[6]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == 0) //low face at z==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[3]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[2]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[3]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == depth-1) //high face at z=depth
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[5]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[6]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
            }//i
        }//j
    }//k
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // three loops, just like the destructor. note we don't delete boundary nodes.
    for (unsigned i=0; i<mBoundaryElements.size(); i++)
    {
        delete mBoundaryElements[i];
    }
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
    mBoundaryElements.clear();
    mDeletedElementIndices.clear();
    mDeletedBoundaryElementIndices.clear();
    mDeletedNodeIndices.clear();
    mBoundaryNodes.clear();
    mAddedNodes = false;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateBoundaryOfFlaggedRegion()
{
    // a set of nodes which lie on the face, size 3 in 2D, size 4 in 3D
    typedef std::set<unsigned> FaceNodes; 

    // face maps to true the first time it is encountered, and false subsequent 
    // times. Thus, faces mapping to true at the end are boundary faces
    std::map<FaceNodes,bool> face_on_boundary;

    // loop over all elements
    ElementIterator iter = GetElementIteratorBegin();
    while (iter != GetElementIteratorEnd())
    {
        if((*iter)->IsFlagged())
        {
            // to get faces, initially start with all nodes..
            std::set<unsigned> all_nodes;
            for(unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                all_nodes.insert( (*iter)->GetNodeGlobalIndex(i) );
            }
            
            // remove one node in turn to obtain each face
            for(unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                FaceNodes face_nodes = all_nodes;
                face_nodes.erase( (*iter)->GetNodeGlobalIndex(i) );

                // search the map of faces to see if it contains this face                 
                std::map<FaceNodes,bool>::iterator it=face_on_boundary.find(face_nodes);

                if(it == face_on_boundary.end())
                {
                    // face not found, add and assume on boundary
                    face_on_boundary[face_nodes]=true;
                }
                else
                {
                    // face found in map, so not on boundary
                    it->second = false;
                }
            }
            
        }
        iter++;
    }
    
    // boundary nodes to be returned
    std::set<unsigned> boundary_of_flagged_region;
    
    // get all faces in the map
    std::map<FaceNodes,bool>::iterator it=face_on_boundary.begin();
    while(it!=face_on_boundary.end())
    {
        // if the face maps to true it is on the boundary
        if(it->second==true)
        {
            // get all nodes in the face and put in set to be returned
            boundary_of_flagged_region.insert(it->first.begin(),it->first.end());
        }
        it++;
    }
        
    return boundary_of_flagged_region;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector = rLocationB - rLocationA;
    
    return vector;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetAngleBetweenNodes(unsigned indexA, unsigned indexB)
{
    assert(SPACE_DIM == 2);
    assert(SPACE_DIM == ELEMENT_DIM);

    double x_diff = mNodes[indexB]->rGetLocation()[0] - mNodes[indexA]->rGetLocation()[0];
    double y_diff = mNodes[indexB]->rGetLocation()[1] - mNodes[indexA]->rGetLocation()[1];
    
    if (x_diff==0)
    {
        if (y_diff>0)
        {
            return M_PI/2.0;
        }
        else if (y_diff<0)
        {
            return -M_PI/2.0;
        }
        else
        {
            EXCEPTION("Tried to compute polar angle of (0,0)");
        }
    } 
    
    double angle = atan2(y_diff,x_diff);
    return angle;    
}


   
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWidth(const unsigned& rDimension) const
{
    assert(rDimension < SPACE_DIM);
    c_vector<double,2> extremes = GetWidthExtremes(rDimension);
    return extremes[1]-extremes[0];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,2> ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWidthExtremes(const unsigned& rDimension) const
{
    assert(rDimension < SPACE_DIM);
    double max = -1e200;
    double min = 1e200;
    assert(GetNumAllNodes() > 0u);
    for (unsigned i=0 ; i<GetNumAllNodes() ; i++)
    {
        if (!mNodes[i]->IsDeleted())
        {
            double this_node_value = mNodes[i]->rGetLocation()[rDimension];
            if (this_node_value>max) 
            {
                max = this_node_value;  
            }
            if (this_node_value < min)
            {
                min = this_node_value;
            }
        }
    }
    c_vector<double,2> extremes;
    extremes[0] = min;
    extremes[1] = max;
    return extremes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::UnflagAllElements()
{
    ElementIterator i_element;
    for (i_element = GetElementIteratorBegin();
         i_element != GetElementIteratorEnd();
         i_element++)
    {
         (*i_element)->Unflag();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::FlagElementsNotContainingNodes(std::set<unsigned> nodesList)
{
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator
        iter = GetElementIteratorBegin();
    while (iter != GetElementIteratorEnd())
    {
        Element<ELEMENT_DIM, SPACE_DIM>& element = **iter;

        bool found_node = false;
        for (unsigned i=0; i<element.GetNumNodes(); i++)
        {
            unsigned node_index = element.GetNodeGlobalIndex(i);
            
            std::set<unsigned>::iterator set_iter = nodesList.find(node_index);
            if(set_iter!=nodesList.end())
            {
                found_node = true;
            }               
        }
        
        if (!found_node)
        {
            element.Flag();
        }
        ++iter;
    }
}


//////////////////////////////////////////////////////////////////////////////
//                          edge iterator class                           // 
//////////////////////////////////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeA()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeALocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeB()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeBLocalIndex);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator!=(const ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& other)
{
    return (mElemIndex != other.mElemIndex ||
            mNodeALocalIndex != other.mNodeALocalIndex ||
            mNodeBLocalIndex != other.mNodeBLocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator++()
{
    std::set<unsigned> current_node_pair;
    std::set<std::set<unsigned> >::iterator set_iter;
    
    do
    {
        // Advance to the next edge in the mesh.
        // Node indices are incremented modulo #nodes_per_elem
        mNodeBLocalIndex = (mNodeBLocalIndex + 1) % (ELEMENT_DIM+1);
        if (mNodeBLocalIndex == mNodeALocalIndex)
        {
            mNodeALocalIndex = (mNodeALocalIndex + 1) % (ELEMENT_DIM+1);
            mNodeBLocalIndex = (mNodeALocalIndex + 1) % (ELEMENT_DIM+1);
        }
        
        if (mNodeALocalIndex == 0 && mNodeBLocalIndex == 1) // advance to next element...
        {
            mElemIndex++;
            // ...skipping deleted ones
            while(mElemIndex!=mrMesh.GetNumAllElements() && mrMesh.GetElement(mElemIndex)->IsDeleted())
            {
                mElemIndex++;
            }
        }

        if(mElemIndex != mrMesh.GetNumAllElements())
        {
            unsigned node_a_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeALocalIndex);
            unsigned node_b_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeBLocalIndex);
        
            // Check we haven't seen it before
            current_node_pair.clear();
            current_node_pair.insert(node_a_global_index);
            current_node_pair.insert(node_b_global_index);
            set_iter = mEdgesVisited.find(current_node_pair);
        } 
    }
    while (*this != mrMesh.EdgesEnd() && set_iter != mEdgesVisited.end());
    mEdgesVisited.insert(current_node_pair);
    
    return (*this);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::EdgeIterator(ConformingTetrahedralMesh& rMesh, unsigned elemIndex)
    : mrMesh(rMesh),
      mElemIndex(elemIndex),
      mNodeALocalIndex(0),
      mNodeBLocalIndex(1)
{
    if(elemIndex==mrMesh.GetNumAllElements())
    {
        return;
    }
    
    mEdgesVisited.clear();
    
    // add the current node pair to the store
    std::set<unsigned> current_node_pair;
    unsigned node_a_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeALocalIndex);
    unsigned node_b_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeBLocalIndex);
    current_node_pair.insert(node_a_global_index);
    current_node_pair.insert(node_b_global_index);
    
    mEdgesVisited.insert(current_node_pair);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesBegin()
{
    unsigned first_element_index=0;
    while(first_element_index!=GetNumAllElements() && GetElement(first_element_index)->IsDeleted())
    {
        first_element_index++;
    }
    return EdgeIterator(*this, first_element_index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesEnd()
{
    return EdgeIterator(*this, GetNumAllElements());
}



#endif //_CONFORMINGTETRAHEDRALMESH_HPP_
