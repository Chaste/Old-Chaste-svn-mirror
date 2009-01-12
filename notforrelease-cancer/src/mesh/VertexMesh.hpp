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
#include <map>
#include <algorithm>

#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>

#include "AbstractMesh.hpp"
#include "VertexMeshReader2d.hpp"
#include "VertexElement.hpp"
#include "BoundaryElement.hpp"
#include "NodeMap.hpp"

/**
 * A vertex-based mesh class, for use in vertex-based tissue simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh
{
    friend class TestVertexMesh;
    
private:

    /** Vector of pointers to nodes. */
    std::vector<Node<SPACE_DIM>*> mNodes;

    /** Vector of pointers to VertexElements. */    
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;
   
    /** The minimum distance apart that two nodes in the mesh can be without causing element rearrangment. */
    double mThresholdDistance;
    
    /** Indices of nodes that have been deleted. These indices can be reused when adding new elements/nodes. */
    std::vector<unsigned> mDeletedNodeIndices;
    
    /** Whether nodes have been added to the mesh. */
    bool mAddedNodes;

    /** Create correspondences between VertexElements and Nodes in the mesh. */
    void SetupVertexElementsOwnedByNodes();
    
    /**
     * Helper method for ReMesh to perform the T1 Swap
     * 
     * @param rNodeA one of the nodes to perform the swap with 
     * @param rNodeB the other node to perform the swap
     */  
    void PerformT1Swap(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mThresholdDistance;      
    }
    
public:

    /**
     * Default constructor.
     * 
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     * @param thresholdDistance the minimum threshold distance for element rearrangment (defaults to 0.01)
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, 
               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
               double thresholdDistance=0.01);

    /**
     * Constructor for use by serializer.
     *
     * @param thresholdDistance the minimum threshold distance for element rearrangment (defaults to 0.01)
     */
    VertexMesh(double thresholdDistance=0.01);
    
    /**
     * Destructor.
     */
    ~VertexMesh();
    
    /**
     * @return mThresholdDistance
     */
    double GetThresholdDistance() const;
    /**
     * Helper constructor, creates a rectangular vertex-based mesh.
     * 
     * @param numAcross number of VertexElements across
     * @param numUp number of VertexElements up
     * @param thresholdDistance the minimum threshold distance for element rearrangment (defaults to 0.01)
     */
    VertexMesh(unsigned numAcross, unsigned numUp, double thresholdDistance=0.01);

    /**
     * @return the number of Nodes in the mesh.
     */
    unsigned GetNumNodes() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    unsigned GetNumElements() const;

    /**
     * @param index the global index of a specified node
     * 
     * @return a pointer to the node
     */
    Node<SPACE_DIM>* GetNode(unsigned index) const;
    
    /**
     * @param index  the global index of a specified vertex element
     * 
     * @return a pointer to the vertex element
     */    
    VertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * Add a node to the mesh.
     *
     * NB. After calling this one or more times, you must then call ReMesh
     *
     */
    unsigned AddNode(Node<SPACE_DIM> *pNewNode);

    /**
     *  Move the node with a particular index to a new point in space.
     * 
      * @param nodeIndex the index of the node to be moved
      * @param point the new target location of the node
      */
    void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point);

    /**
     * Delete mNodes and mElements.
     */
    void Clear();
    
    /**
     * Re-mesh the mesh.
     * 
     * @param map a NodeMap which associates the indices of VertexElements in the old mesh
     *            with indices of VertexElements in the new mesh.  This should be created 
     *            with the correct size, GetNumElements()
     */
    void ReMesh(NodeMap& elementMap);

    /**
     * Alternative version of remesh which takes no parameters does not require a NodeMap. 
     * Note: inherited classes should overload ReMesh(NodeMap&).
     */
    void ReMesh();
       
    /** 
     * Construct the mesh using a mesh reader.
     * 
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(VertexMeshReader2d& rMeshReader);

};

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Cylindrical2dMesh
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const VertexMesh<ELEMENT_DIM, SPACE_DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const double threshold_distance = t->GetThresholdDistance();
    ar << threshold_distance;
}

/**
 * De-serialize constructor parameters and initialise Cylindrical2dMesh.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, VertexMesh<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double threshold_distance;
    ar >> threshold_distance;

    // Invoke inplace constructor to initialise instance
    ::new(t)VertexMesh<ELEMENT_DIM, SPACE_DIM>(threshold_distance);
}
}
} // namespace ...

#endif /*VERTEXMESH_HPP_*/
