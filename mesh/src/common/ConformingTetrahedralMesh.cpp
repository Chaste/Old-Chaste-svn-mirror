#ifndef _CONFORMINGTETRAHEDRALMESH_CPP_
#define _CONFORMINGTETRAHEDRALMESH_CPP_

#include "ConformingTetrahedralMesh.hpp"
#include "Exception.hpp"

#include <vector>
#include <map>
#include <set>
#include <algorithm>

template<int ELEMENT_DIM, int SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConformingTetrahedralMesh()
{}


template<int ELEMENT_DIM, int SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConformingTetrahedralMesh(long numElements)
{
    mElements.reserve(numElements);
}



template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    int orderOfBasisFunctions)
{

    // We only use linear or quadratic basis functions
    assert(orderOfBasisFunctions == 1 || orderOfBasisFunctions == 2);
    
    // Record number of corner nodes
    mNumCornerNodes = rMeshReader.GetNumNodes();
    
    // Reserve memory for nodes, so we don't have problems with pointers stored in
    // elements becoming invalid.
    // For quadratic basis functions we need to go through and count them first
    
    if (orderOfBasisFunctions == 1)
    {
        mNodes.reserve(mNumCornerNodes);
    }
    else if (orderOfBasisFunctions == 2)
    {
        typename std::map<std::pair<int,int>,int>::const_iterator iterator;
        std::map<std::pair<int,int>,int> internal_nodes_map;
        
        std::vector<Node<SPACE_DIM> > temp_nodes;
        temp_nodes.reserve(mNumCornerNodes);
        
        // Add corner nodes to temporary nodes vector
        std::vector<double> coords;
        for (int i=0; i < mNumCornerNodes; i++)
        {
            coords = rMeshReader.GetNextNode();
            temp_nodes.push_back(Node<SPACE_DIM>(i, Point<SPACE_DIM>(coords), false));
        }
        
        int new_node_index = mNumCornerNodes;
        
        // Create the map of internal node indices
        for (int i=0; i < rMeshReader.GetNumElements(); i++)
        {
            std::vector<int> node_indices = rMeshReader.GetNextElement();
            std::vector<const Node<SPACE_DIM>*> nodes;
            
            //Debugging code
            int nodes_size = node_indices.size();
            
            for (int j=0; j < nodes_size; j++)
            {
                assert(node_indices[j] < (int) temp_nodes.size());
                nodes.push_back(&temp_nodes[node_indices[j]]);
            }
            
            for (int j=0; j < ELEMENT_DIM + 1; j++)
            {
                for (int k=j+1; k < ELEMENT_DIM + 1; k++)
                {
                    assert(j < nodes_size && k < nodes_size);
                    int node_i = nodes[j]->GetIndex();
                    int node_j = nodes[k]->GetIndex();
                    if (node_j < node_i)
                    {
                        int temp = node_i;
                        node_i = node_j;
                        node_j = temp;
                    }
                    iterator = internal_nodes_map.find(std::pair<int,int>(node_i, node_j));
                    if (iterator == internal_nodes_map.end())
                    {
                        // add node to map
                        internal_nodes_map[(std::pair<int,int>(node_i, node_j))] = new_node_index;
                        new_node_index++;
                    }
                }
            }
        }
        mNodes.reserve(new_node_index);
    }
    
    rMeshReader.Reset();
    
    typename std::map<std::pair<int,int>,int>::const_iterator iterator;
    std::map<std::pair<int,int>,int> internal_nodes_map;
    
    // Add corner nodes
    std::vector<double> coords;
    for (int i=0; i < mNumCornerNodes; i++)
    {
        coords = rMeshReader.GetNextNode();
        mNodes.push_back(new Node<SPACE_DIM>(i, Point<SPACE_DIM>(coords), false));
    }
    
    int new_node_index = mNumCornerNodes;
    
    if (orderOfBasisFunctions == 2)
    {
        for (int i=0; i < rMeshReader.GetNumElements(); i++)
        {
            std::vector<int> node_indices = rMeshReader.GetNextElement();
            std::vector<const Node<SPACE_DIM>*> nodes;
            
            //Debugging code
            int nodes_size = node_indices.size();
            
            for (int j=0; j < nodes_size; j++)
            {
                assert(node_indices[j] < (int) mNodes.size());
                nodes.push_back(mNodes[node_indices[j]]);
            }
            
            for (int j=0; j < ELEMENT_DIM + 1; j++)
            {
                for (int k=j+1; k < ELEMENT_DIM + 1; k++)
                {
                    assert(j < nodes_size && k < nodes_size);
                    int node_i = nodes[j]->GetIndex();
                    int node_j = nodes[k]->GetIndex();
                    if (node_j < node_i)
                    {
                        int temp = node_i;
                        node_i = node_j;
                        node_j = temp;
                    }
                    iterator = internal_nodes_map.find(std::pair<int,int>(node_i, node_j));
                    if (iterator == internal_nodes_map.end())
                    {
                        // add node to map
                        internal_nodes_map[(std::pair<int,int>(node_i, node_j))] = new_node_index;
                        // add node to mesh
                        const Node<SPACE_DIM>* node1 = GetNodeAt(node_i);
                        const Node<SPACE_DIM>* node2 = GetNodeAt(node_j);
                        Node<SPACE_DIM> *p_new_node=new Node<SPACE_DIM>(new_node_index,
                                                                        node1->GetPoint().MidPoint(node2->GetPoint()),
                                                                        node1->IsBoundaryNode() && node2->IsBoundaryNode());
                        mNodes.push_back(p_new_node);
                        new_node_index++;
                    }
                }
            }
        }
    }
    rMeshReader.Reset();
    // Add elements
    new_node_index = mNumCornerNodes;
    mElements.reserve(rMeshReader.GetNumElements());
    
    for (unsigned element_index=0; element_index < (unsigned) rMeshReader.GetNumElements(); element_index++)
    {
        std::vector<int> node_indices = rMeshReader.GetNextElement();
        std::vector<Node<SPACE_DIM>*> nodes;
        int nodes_size = node_indices.size();
        for (int j=0; j<nodes_size; j++)
        {
            assert(node_indices[j] < (int) mNodes.size());
            nodes.push_back(mNodes[node_indices[j]]);
        }
        
        if (orderOfBasisFunctions == 2)
        {
            for (int j=0; j < ELEMENT_DIM + 1; j++)
            {
                for (int k=j+1; k < ELEMENT_DIM + 1; k++)
                {
                    assert(j < nodes_size && k < nodes_size);
                    int node_i = nodes[j]->GetIndex();
                    int node_j = nodes[k]->GetIndex();
                    if (node_j < node_i)
                    {
                        int temp = node_i;
                        node_i = node_j;
                        node_j = temp;
                    }
                    iterator = internal_nodes_map.find(std::pair<int,int>(node_i, node_j));
                    assert(iterator != internal_nodes_map.end());
                    // add node to element
                    nodes.push_back(this->GetNodeAt(iterator->second));
                    new_node_index++;
                }
            }
        }
        mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes,orderOfBasisFunctions));
    }
    
    // Add boundary elements & nodes
    int actual_face_index=0;
    for (unsigned face_index=0; face_index<(unsigned)rMeshReader.GetNumFaces(); face_index++)
    {
        std::vector<int> node_indices = rMeshReader.GetNextFace();
        
        // Determine if this is a boundary face
        std::set<unsigned> containing_element_indices; // Elements that contain this face
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned node_index=0; node_index<node_indices.size(); node_index++)
        {
            assert(node_indices[node_index] < (int) mNodes.size());
            // Add Node pointer to list for creating an element
            nodes.push_back(mNodes[node_indices[node_index]]);
            
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
        
        //If the following assertion is thrown, it means that the .edge/.face file does not
        //match the .ele file -- they were generated at separate times.  Simply remove the internal
        //edges/faces by hand.
        assert(containing_element_indices.size() != 0);
        
        if (containing_element_indices.size() == 1)
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
            
            
            if (orderOfBasisFunctions == 2)
            {
                int nodes_size = node_indices.size();
                for (int j=0; j < ELEMENT_DIM; j++)
                {
                    for (int k=j+1; k < ELEMENT_DIM; k++)
                    {
                        assert(j < nodes_size && k < nodes_size);
                        int node_i = nodes[j]->GetIndex();
                        int node_j = nodes[k]->GetIndex();
                        if (node_j < node_i)
                        {
                            int temp = node_i;
                            node_i = node_j;
                            node_j = temp;
                        }
                        iterator = internal_nodes_map.find(std::pair<int,int>(node_i, node_j));
                        assert(iterator != internal_nodes_map.end());
                        // add node to element
                        nodes.push_back(this->GetNodeAt(iterator->second));
                    }
                }
            }
            
            // The added elements will be deleted in our destructor
            mBoundaryElements.push_back(
                new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(actual_face_index,nodes,orderOfBasisFunctions));
            actual_face_index++;
            
        }
    }
}

template<int ELEMENT_DIM, int SPACE_DIM>
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


template<int ELEMENT_DIM, int SPACE_DIM>
int ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM> *pNewNode)
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
    
    return pNewNode->GetIndex();
}


/**
 * Get a node reference from the mesh.
 *
 * Note that this may become invalid if nodes are subsequently added to the mesh.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
Node<SPACE_DIM> *ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNodeAt(long index)
{
    assert(index < (long) mNodes.size());
    return (mNodes[index]);
}

/// Returns the number of nodes that are actually in use
template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNodes.size() - mDeletedNodeIndices.size();
}
template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes()
{
    return mNodes.size();
}

template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mElements.size() - mDeletedElementIndices.size();
}
template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements()
{
    return mElements.size();
}

template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryNodes()
{
    return mBoundaryNodes.size();
}

template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements()
{
    return mBoundaryElements.size() - mDeletedBoundaryElementIndices.size();
}

template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllBoundaryElements()
{
    return mBoundaryElements.size();
}


template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumCornerNodes()
{
    return mNumCornerNodes;
}



template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RescaleMeshFromBoundaryNode(Point<1> updatedPoint, int boundaryNodeIndex)
{

    double scaleFactor = updatedPoint[0] / GetNodeAt(boundaryNodeIndex)->GetPoint()[0];
    double temp;
    for (int i=0; i < boundaryNodeIndex+1; i++)
    {
        temp = scaleFactor * mNodes[i]->GetPoint()[0];
        Point<1> newPoint(temp);
        mNodes[i]->SetPoint(newPoint);
    }
    
}

template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned index, Point<SPACE_DIM> point, bool verify)
{
    mNodes[index]->SetPoint(point);
    if (verify)
    {
        for (int i=0; i<mNodes[index]->GetNumContainingElements(); i++)
        {
            try
            {
                GetElement(mNodes[index]->GetNextContainingElementIndex())->RefreshJacobianDeterminant();
            }
            catch (Exception e)
            {
                EXCEPTION("Moving node caused an element to have a non-positive Jacobian determinant");
            }
        }
        for (int i=0; i<mNodes[index]->GetNumBoundaryElements(); i++)
        {
            try
            {
                GetBoundaryElement(mNodes[index]->GetNextBoundaryElementIndex())->RefreshJacobianDeterminant();
            }
            catch (Exception e)
            {
                EXCEPTION("Moving node caused a boundary element to have a non-positive Jacobian determinant");
            }
        }
    }
}

/**
 * This method allows the mesh properties to be re-calculated after one
 * or more node have been moved.
 *
 */
template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
    for (unsigned i=0; i<mElements.size();i++)
    {
        mElements[i]->RefreshJacobianDeterminant();
    }
    
    //Refresh each boundary element
    for (unsigned i=0; i<mBoundaryElements.size();i++)
    {
        mBoundaryElements[i]->RefreshJacobianDeterminant();
    }
    
}

template<int ELEMENT_DIM, int SPACE_DIM>
int ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefineElement(
    Element<ELEMENT_DIM,SPACE_DIM>* pElement,
    Point<SPACE_DIM> point)
{

    // Add a new node from the point that is passed to RefineElement
    int new_node_index = AddNode(new Node<SPACE_DIM>(0, point));
    // Note: the first argument is the index of the node, which is going to be
    //       overriden by AddNode, so it can safely be ignored
    
    //This loop constructs the extra elements which are going to fill the space
    for (unsigned i = 0; i < ELEMENT_DIM; i++)
    {
        // First, make a copy of the current element making sure we update its index
        int new_elt_index;
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
        
        try
        {
            p_new_element->RefreshJacobianDeterminant();
        }
        catch (Exception e)
        {
            //Clean up this element, because it's not going to be added to the mesh
            delete p_new_element;
            EXCEPTION("RefineElement could not be completed (point was not in element)");
        }
        
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
    try
    {
        pElement->RefreshJacobianDeterminant();
    }
    catch (Exception e)
    {
        //There's no memory management to do here, because this original element is
        //already referenced by the mesh
        EXCEPTION("RefineElement could not be completed (point was not in element)");
    }
    
    
    return new_node_index;
}


template<int ELEMENT_DIM, int SPACE_DIM>
double ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMeshVolume()
{
    double mesh_volume = 0.0;
    
    ElementIterator it = GetElementIteratorBegin();
    
    while (it != GetElementIteratorEnd())
    {
        mesh_volume += (*it)->GetJacobianDeterminant();
        it++;
    }
    
    if (ELEMENT_DIM == 2)
    {
        mesh_volume /= 2.0;
    }
    else if (ELEMENT_DIM == 3)
    {
        mesh_volume /= 6.0;
    }
    return mesh_volume;
}

template<int ELEMENT_DIM, int SPACE_DIM>
double ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMeshSurface()
{
    //ELEMENT_DIM-1 is the dimension of the boundary element
    int bound_element_dim=ELEMENT_DIM-1;
    assert(bound_element_dim < 3 &&  bound_element_dim>=0);
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
 * Translate the mesh.
 * @param xMovement is the x-displacement,
 * @param yMovement is the y-displacement,
 * @param zMovement is the z-displacement,
 **/
template <int ELEMENT_DIM, int SPACE_DIM>
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
template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Translate(c_vector<double, SPACE_DIM> transVec)
{
    long num_nodes=GetNumAllNodes();
    
    for (int i=0; i<num_nodes; i++)
    {
        Point<SPACE_DIM> point = mNodes[i]->rGetPoint();
        point.rGetLocation() += transVec;
        mNodes[i]->SetPoint(point);
    }
    
    RefreshMesh();
}




/**
 * Do a general mesh rotation with an +ve determinant orthonormal rotation_matrix
 * This is the method that actually does the work
 * @param rotation_matrix is a Ublas rotation matrix of the correct form
 **/
template <int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(
    c_matrix<double , SPACE_DIM, SPACE_DIM> rotation_matrix)
{
    long num_nodes=GetNumAllNodes();
    for (int i=0; i<num_nodes; i++)
    {
        Point<SPACE_DIM> point = mNodes[i]->rGetPoint();
        point.rGetLocation() = prod(rotation_matrix, point.rGetLocation());
        mNodes[i]->SetPoint(point);
    }
    
    RefreshMesh();
}

/**
* Do an angle axis rotation
* @param axis is the axis of rotation (does not need to be normalised)
* @param angle is the angle in radians
 **/
template<int ELEMENT_DIM, int SPACE_DIM>
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
template <int ELEMENT_DIM, int SPACE_DIM>
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
template <int ELEMENT_DIM, int SPACE_DIM>
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
template <int ELEMENT_DIM, int SPACE_DIM>
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


template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DeleteBoundaryNodeAt(long index)
{
    if (!mNodes[index]->IsBoundaryNode() )
    {
        EXCEPTION(" You may only delete a boundary node ");
    }
    
    mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);
    // Update the boundary node vector
    typename std::vector<const Node<SPACE_DIM>*>::iterator b_node_iter
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
        for (int i=0 ; i< p_element->GetNumNodes();i++)
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

#endif // _CONFORMINGTETRAHEDRALMESH_CPP_
