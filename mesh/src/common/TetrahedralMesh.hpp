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


#ifndef _TETRAHEDRALMESH_HPP_
#define _TETRAHEDRALMESH_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>

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
class TetrahedralMesh
{
    friend class TestTetrahedralMesh; // to give access to private methods (not variables)
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

    std::vector< Node<SPACE_DIM> *> mBoundaryNodes;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       // Don't do anything - this is just so subclasses can archive member variables.
    }        

public:

    TetrahedralMesh();
    TetrahedralMesh(unsigned numElements);
    //TetrahedralMesh(std::vector<Node<SPACE_DIM> *> nodes);

    virtual ~TetrahedralMesh();

    virtual void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                 bool cullInternalFaces=false);

    Node<SPACE_DIM> *GetNode(unsigned index);

    virtual unsigned GetNumNodes();
    virtual unsigned GetNumElements();
    virtual unsigned GetNumBoundaryElements();
    unsigned GetNumAllNodes() const;
    unsigned GetNumAllElements();
    unsigned GetNumAllBoundaryElements();
    unsigned GetNumBoundaryNodes();// should this be overloaded and virtual too?

    void ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile);
    std::vector<unsigned>& rGetNodesPerProcessor();

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

    void RefreshMesh(void);

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
    virtual void Clear();

    /**
     *  Return the set of nodes which are on the boundary of the flagged region(s)
     */
    std::set<unsigned> CalculateBoundaryOfFlaggedRegion();

    /**
     * Returns distance between two nodes
     *
     * @param indexA a node index
     * @param indexB a node index
     *
     * @return straight line distance between two nodes.
     *
     * N.B. This calls GetDistanceBetweenNodes which can be overridden 
     * in daughter classes e.g. Cylindrical2dMesh.  Therefore the distance
     * is not necessarily Euclidean
     *
     */
    double GetDistanceBetweenNodes(unsigned indexA, unsigned indexB);

    /**
     * Returns a vector between two points in space
     *
     * @param rLocationA a c_vector of co-ordinates
     * @param rLocationB a c_vector of co-ordinates
     *
     * @return vector from location A to location B.
     *
     * N.B. This can be overridden in daughter classes
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
        EdgeIterator(TetrahedralMesh& rMesh, unsigned elemIndex);

    private:
        /** Keep track of what edges have been visited */
        std::set<std::set<unsigned> > mEdgesVisited;

        TetrahedralMesh& mrMesh;

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
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::TetrahedralMesh()
{
    Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::TetrahedralMesh(unsigned numElements)
{
    Clear();
    mElements.reserve(numElements);
}

//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::TetrahedralMesh(std::vector<Node<SPACE_DIM> *> nodes)
//  //: mNodes(nodes)
//{
//    Clear();
//    for (unsigned index=0; index<nodes.size(); index++)
//    {
//        Node<SPACE_DIM>* temp_node = nodes[index];
//        mNodes.push_back(temp_node);
//    }
//    mAddedNodes = true;
//    NodeMap node_map(nodes.size());
//    ReMesh(node_map);
//}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
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

// NOTE: currently just reading element vertices from mesh reader - even if it
// does contain information about internal nodes (ie for quadratics) this is
// ignored here and used elsewhere

        // unsigned nodes_size = node_indices.size();
        for (unsigned j=0; j<ELEMENT_DIM+1; j++) // num vertices=ELEMENT_DIM+1, may not be equal to nodes_size.
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
            // only if not 1D as this assertion does not apply to quadratic 1D meshes
            if(ELEMENT_DIM!=1)
            {
                //If the following assertion is thrown, it means that the .edge/.face file does not
                //match the .ele file -- they were generated at separate times.  Simply remove the internal
                //edges/faces by hand.
                assert(containing_element_indices.size() != 0);
            }

            // if num_containing_elements is greater than 1, it is not an boundary face
            if(containing_element_indices.size() > 1)
            {
                is_boundary_face = false;
            }
            
            // in 1D QUADRATICS, all nodes are faces, so internal nodes which don't have any
            // containing elements must also be unmarked as a boundary face
            if( (ELEMENT_DIM==1) && (containing_element_indices.size()==0))
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
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~TetrahedralMesh()
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile)
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
std::vector<unsigned>& TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodesPerProcessor()
{
    return mNodesPerProcessor;
}

/**
 * Get a node reference from the mesh.
 *
 * Note that this may become invalid if nodes are subsequently added to the mesh.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index)
{
    assert(index < mNodes.size());
    return (mNodes[index]);
}

/// Returns the number of nodes that are actually in use
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes() const
{
    return mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements()
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryNodes()
{
    return mBoundaryNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements()
{
    return mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllBoundaryElements()
{
    return mBoundaryElements.size();
}








template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateVolume()
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
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateSurfaceArea()
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Scale(
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Translate(
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Translate(c_vector<double, SPACE_DIM> transVec)
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(c_vector<double,3> axis, double angle)
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateX(const double theta)
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateY(const double theta)
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateZ(const double theta)
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



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes(std::vector<unsigned>& perm)
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodesWithMetisBinaries(unsigned numProcs)
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

    /*
     *  Create the permutation vector based on Metis output
     */
    std::vector<unsigned> permutation(GetNumNodes(), UINT_MAX);
    std::vector<unsigned> count(numProcs,0u);

    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        unsigned part = partition[node_index];
        // Permutation defined like: new index for node node_index is "offset[part] + count[part]"
        permutation [ node_index ] = offset[part] + count[part];

        count[part]++;
    }

    PermuteNodes(permutation);

}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructLinearMesh(unsigned width)
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRectangularMesh(unsigned width, unsigned height, bool stagger)
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
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndex(ChastePoint<SPACE_DIM> testPoint, bool strict, std::set<unsigned> testElements)
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
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNearestElementIndex(ChastePoint<SPACE_DIM> testPoint)
{
    ///\todo This ought to return a set of all elements that contain the point (if the point is a node in the mesh then it's contained in multiple elements)
    ///\todo Polling every element is unnecessary.  We ought to start from a likely place and hill climb

    double max_min_weight=-INFINITY;
    unsigned closest_index=0;
    for (unsigned i=0; i < mElements.size();i++)
    {
        ///\todo What if the element is deleted?
        c_vector<double, ELEMENT_DIM+1> weight=mElements[i]->CalculateInterpolationWeights(testPoint);
        double neg_weight_sum=0.0;
        for (unsigned j=0; j<=ELEMENT_DIM; j++)
        {
            if (weight[j]<0.0)
            {
                neg_weight_sum+=weight[j];
            }
        }
        if (neg_weight_sum > max_min_weight)
        {
            max_min_weight = neg_weight_sum;
            closest_index=i;
        }

    }
    return closest_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndices(ChastePoint<SPACE_DIM> testPoint)
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructCuboid(unsigned width,
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
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
    mBoundaryNodes.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateBoundaryOfFlaggedRegion()
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
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetDistanceBetweenNodes(unsigned indexA, unsigned indexB)
{
    c_vector<double, SPACE_DIM> vector = GetVectorFromAtoB(mNodes[indexA]->rGetLocation(), mNodes[indexB]->rGetLocation());
    return norm_2(vector);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector = rLocationB - rLocationA;

    return vector;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetAngleBetweenNodes(unsigned indexA, unsigned indexB)
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
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWidth(const unsigned& rDimension) const
{
    assert(rDimension < SPACE_DIM);
    c_vector<double,2> extremes = GetWidthExtremes(rDimension);
    return extremes[1]-extremes[0];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,2> TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWidthExtremes(const unsigned& rDimension) const
{
    assert(rDimension < SPACE_DIM);
    double max = -1e200;
    double min = 1e200;
    assert(GetNumAllNodes() > 0u);
    for (unsigned i=0; i<GetNumAllNodes(); i++)
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::UnflagAllElements()
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
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::FlagElementsNotContainingNodes(std::set<unsigned> nodesList)
{
    typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator
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
Node<SPACE_DIM>* TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeA()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeALocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeB()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeBLocalIndex);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator!=(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& other)
{
    return (mElemIndex != other.mElemIndex ||
            mNodeALocalIndex != other.mNodeALocalIndex ||
            mNodeBLocalIndex != other.mNodeBLocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator++()
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
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::EdgeIterator(TetrahedralMesh& rMesh, unsigned elemIndex)
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
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesBegin()
{
    unsigned first_element_index=0;
    while(first_element_index!=GetNumAllElements() && GetElement(first_element_index)->IsDeleted())
    {
        first_element_index++;
    }
    return EdgeIterator(*this, first_element_index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesEnd()
{
    return EdgeIterator(*this, GetNumAllElements());
}

/**
 * This method allows the mesh properties to be re-calculated after one
 * or more node have been moved.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
    for (unsigned i=0; i<mElements.size();i++)
    {
        if (!mElements[i]->IsDeleted())
        {
            mElements[i]->RefreshJacobianDeterminant();
        }
    }

    //Refresh each boundary element
    for (unsigned i=0; i<this->mBoundaryElements.size();i++)
    {
        if (!this->mBoundaryElements[i]->IsDeleted())
        {
            try
            {
                this->mBoundaryElements[i]->RefreshJacobianDeterminant();
            }
            catch (Exception e)
            {
                //Since we may have rotated the mesh, it's okay for normals to swing round
            }
        }
    }

}



#endif //_TETRAHEDRALMESH_HPP_
