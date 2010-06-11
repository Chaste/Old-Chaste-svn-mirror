/*

Copyright (C) University of Oxford, 2005-2010

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

#include "VertexMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include <list>

/**
 * Global method allowing alist of pairs (unsigned, double) to be compared
 * in terms of their second entry and std::list.sort() to be called.
 */
bool IndexAngleComparison(const std::pair<unsigned, double> lhs, const std::pair<unsigned, double> rhs)
{
    return lhs.second < rhs.second;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                               std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertexElements)
    : mpDelaunayMesh(NULL)
{

    // Reset member variables and clear mNodes and mElements
    Clear();

    // Populate mNodes and mElements
    for (unsigned node_index=0; node_index<nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index=0; elem_index<vertexElements.size(); elem_index++)
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = vertexElements[elem_index];
        mElements.push_back(p_temp_vertex_element);
    }

    // In 3D, populate mFaces
    if (SPACE_DIM == 3)
    {
        // Use a std::set to keep track of which faces have been added to mFaces
        std::set<unsigned> faces_counted;

        // Loop over mElements
        for (unsigned elem_index=0; elem_index<mElements.size(); elem_index++)
        {
            // Loop over faces of this element
            for (unsigned face_index=0; face_index<mElements[elem_index]->GetNumFaces(); face_index++)
            {
                VertexElement<ELEMENT_DIM-1, SPACE_DIM>* p_face = mElements[elem_index]->GetFace(face_index);

                // If this face is not already contained in mFaces, add it, and update faces_counted
                if (faces_counted.find(p_face->GetIndex()) == faces_counted.end())
                {
                    mFaces.push_back(p_face);
                    faces_counted.insert(p_face->GetIndex());
                }
            }
        }
    }

    // Register elements with nodes
    for (unsigned index=0; index<mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = mElements[index];
        for (unsigned node_index=0; node_index<p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->mMeshChangesDuringSimulation = false;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                           std::vector<VertexElement<ELEMENT_DIM-1, SPACE_DIM>*> faces,
                           std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements)
    : mpDelaunayMesh(NULL)
{
    // Reset member variables and clear mNodes, mFaces and mElements
    Clear();

    // Populate mNodes mFaces and mElements
    for (unsigned node_index=0; node_index<nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }

    for (unsigned face_index=0; face_index<faces.size(); face_index++)
    {
        VertexElement<ELEMENT_DIM-1, SPACE_DIM>* p_temp_face = faces[face_index];
        mFaces.push_back(p_temp_face);
    }

    for (unsigned elem_index=0; elem_index<vertexElements.size(); elem_index++)
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = vertexElements[elem_index];
        mElements.push_back(p_temp_vertex_element);
    }

    // Register elements with nodes
    for (unsigned index=0; index<mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = mElements[index];
        for (unsigned node_index=0; node_index<p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->mMeshChangesDuringSimulation = false;
}


/**
 * This VertexMesh constructor is currently only defined for 2D meshes.
 *
 * @param rMesh a tetrahedral mesh
 */
template<>
VertexMesh<2,2>::VertexMesh(TetrahedralMesh<2,2>& rMesh)
    : mpDelaunayMesh(&rMesh)
{
    // Reset member variables and clear mNodes, mFaces and mElements
    Clear();

    unsigned num_elements = mpDelaunayMesh->GetNumAllNodes();
    unsigned num_nodes = mpDelaunayMesh->GetNumAllElements();

    // Allocate memory for mNodes and mElements
    this->mNodes.reserve(num_nodes);

    // Create as many elements as there are nodes in the mesh
    mElements.reserve(num_elements);
    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_index);
        mElements.push_back(p_element);
    }

    // Populate mNodes
    GenerateVerticesFromElementCircumcentres(rMesh);

    // Loop over elements of the Delaunay mesh
    for (unsigned i=0; i<num_nodes; i++)
    {
        // Loop over nodes owned by this element in the Delaunay mesh
        for (unsigned local_index=0; local_index<3; local_index++)
        {
            unsigned elem_index = mpDelaunayMesh->GetElement(i)->GetNodeGlobalIndex(local_index);
            unsigned num_nodes_in_elem = mElements[elem_index]->GetNumNodes();
            unsigned end_index = num_nodes_in_elem>0 ? num_nodes_in_elem-1 : 0;

            mElements[elem_index]->AddNode(end_index, this->mNodes[i]);
        }
    }

    // Reorder mNodes anticlockwise
    for (unsigned elem_index=0; elem_index<mElements.size(); elem_index++)
    {
        /**
         * Create a std::list of pairs, where each pair comprises the angle
         * between the centre of the Voronoi element and each node with that
         * node's global index in the Voronoi mesh.
         */
        std::list<std::pair<unsigned, double> > index_angle_list;
        for (unsigned local_index=0; local_index<mElements[elem_index]->GetNumNodes(); local_index++)
        {
            c_vector<double, 2> centre_to_vertex = mpDelaunayMesh->GetVectorFromAtoB(mpDelaunayMesh->GetNode(elem_index)->rGetLocation(),
                                                                                     mElements[elem_index]->GetNodeLocation(local_index));

            double angle = atan2(centre_to_vertex(1), centre_to_vertex(0));
            unsigned global_index = mElements[elem_index]->GetNodeGlobalIndex(local_index);

            std::pair<unsigned, double> pair(global_index, angle);
            index_angle_list.push_back(pair);
        }

        // Sort the list in order of increasing angle
        index_angle_list.sort(IndexAngleComparison);

        // Create a new Voronoi element and pass in the appropriate Nodes, ordered anticlockwise
        VertexElement<2,2>* p_new_element = new VertexElement<2,2>(elem_index);
        unsigned count = 0;
        for (std::list<std::pair<unsigned, double> >::iterator list_iter = index_angle_list.begin();
             list_iter != index_angle_list.end();
             ++list_iter)
        {
            unsigned local_index = count>1 ? count-1 : 0;
            p_new_element->AddNode(local_index, mNodes[list_iter->first]);
            count++;
        }

        // Replace the relevant member of mElements with this Voronoi element
        delete mElements[elem_index];
        mElements[elem_index] = p_new_element;
    }

    this->mMeshChangesDuringSimulation = false;
}


/**
 * This VertexMesh constructor is currently only defined for 3D meshes.
 *
 * @param rMesh a tetrahedral mesh
 */
template<>
VertexMesh<3,3>::VertexMesh(TetrahedralMesh<3,3>& rMesh)
    : mpDelaunayMesh(&rMesh)
{
    // Reset member variables and clear mNodes, mFaces and mElements
    Clear();

    unsigned num_nodes = mpDelaunayMesh->GetNumAllElements();

    // Allocate memory for mNodes
    this->mNodes.reserve(num_nodes);

    // Populate mNodes
    GenerateVerticesFromElementCircumcentres(rMesh);

    std::map<unsigned, VertexElement<3,3>*> index_element_map;
    unsigned face_index = 0;
    unsigned element_index = 0;

    // Loop over each edge of the Delaunay mesh and populate mFaces and mElements
    for (TetrahedralMesh<3,3>::EdgeIterator edge_iterator = mpDelaunayMesh->EdgesBegin();
         edge_iterator != mpDelaunayMesh->EdgesEnd();
         ++edge_iterator)
    {
        Node<3>* p_node_a = edge_iterator.GetNodeA();
        Node<3>* p_node_b = edge_iterator.GetNodeB();

        if ( !(p_node_a->IsBoundaryNode() && p_node_b->IsBoundaryNode()) )
        {
            std::set<unsigned>& node_a_element_indices = p_node_a->rGetContainingElementIndices();
            std::set<unsigned>& node_b_element_indices = p_node_b->rGetContainingElementIndices();
            std::set<unsigned> edge_element_indices;

            std::set_intersection(node_a_element_indices.begin(),
                                  node_a_element_indices.end(),
                                  node_b_element_indices.begin(),
                                  node_b_element_indices.end(),
                                  std::inserter(edge_element_indices, edge_element_indices.begin()));

            c_vector<double,3> edge_vector = p_node_b->rGetLocation() - p_node_a->rGetLocation();
            c_vector<double,3> mid_edge = edge_vector*0.5 + p_node_a->rGetLocation();

            unsigned element0_index = *(edge_element_indices.begin());

            c_vector<double,3> basis_vector1 = mNodes[element0_index]->rGetLocation() - mid_edge;

            c_vector<double,3> basis_vector2;
            basis_vector2[0] = edge_vector[1]*basis_vector1[2] - edge_vector[2]*basis_vector1[1];
            basis_vector2[1] = edge_vector[2]*basis_vector1[0] - edge_vector[0]*basis_vector1[2];
            basis_vector2[2] = edge_vector[0]*basis_vector1[1] - edge_vector[1]*basis_vector1[0];

            /**
             * Create a std::list of pairs, where each pair comprises the angle
             * between the centre of the Voronoi element and each node with that
             * node's global index in the Voronoi mesh.
             */
            std::list<std::pair<unsigned, double> > index_angle_list;

            // Loop over each element containing this edge (i.e. those containing both nodes of the edge)
            for (std::set<unsigned>::iterator index_iter = edge_element_indices.begin();
                 index_iter != edge_element_indices.end();
                 ++index_iter)
            {
                // Calculate angle
                c_vector<double, 3> vertex_vector =  mNodes[*index_iter]->rGetLocation() - mid_edge;

                double local_vertex_dot_basis_vector1 = inner_prod(vertex_vector, basis_vector1);
                double local_vertex_dot_basis_vector2 = inner_prod(vertex_vector, basis_vector2);

                double angle = atan2(local_vertex_dot_basis_vector2, local_vertex_dot_basis_vector1);

                std::pair<unsigned, double> pair(*index_iter, angle);
                index_angle_list.push_back(pair);
            }

            // Sort the list in order of increasing angle
            index_angle_list.sort(IndexAngleComparison);

            // Create face
            VertexElement<2,3>* p_face = new VertexElement<2,3>(face_index);
            face_index++;
            unsigned count = 0;
            for (std::list<std::pair<unsigned, double> >::iterator list_iter = index_angle_list.begin();
                 list_iter != index_angle_list.end();
                 ++list_iter)
            {
                unsigned local_index = count>1 ? count-1 : 0;
                p_face->AddNode(local_index, mNodes[list_iter->first]);
                count++;
            }

            // Add face to list of faces
            mFaces.push_back(p_face);

            // Add face to appropriate elements
            if (!p_node_a->IsBoundaryNode())
            {
                if (index_element_map[p_node_a->GetIndex()])
                {
                    // If there is already an element with this index, add the face to it...
                    index_element_map[p_node_a->GetIndex()]->AddFace(p_face);
                }
                else
                {
                    // ...otherwise create an element, add the face to it, and add to the map
                    mVoronoiElementIndexMap[p_node_a->GetIndex()] = element_index;
                    VertexElement<3,3>* p_element = new VertexElement<3,3>(element_index);
                    element_index++;
                    p_element->AddFace(p_face);
                    index_element_map[p_node_a->GetIndex()] = p_element;
                }
            }
            if (!p_node_b->IsBoundaryNode())
            {
                if (index_element_map[p_node_b->GetIndex()])
                {
                    // If there is already an element with this index, add the face to it...
                    index_element_map[p_node_b->GetIndex()]->AddFace(p_face);
                }
                else
                {
                    // ...otherwise create an element, add the face to it, and add to the map
                    mVoronoiElementIndexMap[p_node_b->GetIndex()] = element_index;
                    VertexElement<3,3>* p_element = new VertexElement<3,3>(element_index);
                    element_index++;
                    p_element->AddFace(p_face);
                    index_element_map[p_node_b->GetIndex()] = p_element;
                }
            }
        }
    }

    // Populate mElements
    unsigned elem_count = 0;
    for (std::map<unsigned, VertexElement<3,3>*>::iterator element_iter = index_element_map.begin();
         element_iter != index_element_map.end();
         ++element_iter)
    {
        mElements.push_back(element_iter->second);
        mVoronoiElementIndexMap[element_iter->first] = elem_count;
        elem_count++;
    }

    this->mMeshChangesDuringSimulation = false;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::GenerateVerticesFromElementCircumcentres(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;
    c_matrix<double, ELEMENT_DIM, SPACE_DIM> inverse_jacobian;
    double jacobian_det;

    // Loop over elements of the Delaunay mesh and populate mNodes
    for (unsigned i=0; i<rMesh.GetNumElements(); i++)
    {
        // Calculate the circumcentre of this element in the Delaunay mesh
        rMesh.GetInverseJacobianForElement(i, jacobian, jacobian_det, inverse_jacobian);
        c_vector<double, SPACE_DIM+1> circumsphere = rMesh.GetElement(i)->CalculateCircumsphere(jacobian, inverse_jacobian);

        c_vector<double, SPACE_DIM> circumcentre;
        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            circumcentre(j) = circumsphere(j);
        }

        // Create a node in the Voronoi mesh at the location of this circumcentre
        this->mNodes.push_back(new Node<SPACE_DIM>(i, circumcentre));
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetEdgeLength(unsigned elementIndex1, unsigned elementIndex2)
{
    assert(SPACE_DIM==2);

    std::set<unsigned> node_indices_1;
    for (unsigned i=0; i<mElements[elementIndex1]->GetNumNodes(); i++)
    {
        node_indices_1.insert(mElements[elementIndex1]->GetNodeGlobalIndex(i));
    }
    std::set<unsigned> node_indices_2;
    for (unsigned i=0; i<mElements[elementIndex2]->GetNumNodes(); i++)
    {
        node_indices_2.insert(mElements[elementIndex2]->GetNodeGlobalIndex(i));
    }

    std::set<unsigned> shared_nodes;
    std::set_intersection(node_indices_1.begin(), node_indices_1.end(),
                          node_indices_2.begin(), node_indices_2.end(),
                          std::inserter(shared_nodes, shared_nodes.begin()));

    assert(shared_nodes.size() == 2);

    unsigned index1 = *(shared_nodes.begin());
    unsigned index2 = *(++(shared_nodes.begin()));

    c_vector<double, SPACE_DIM> node1_location = this->mNodes[index1]->rGetLocation();
    c_vector<double, SPACE_DIM> node2_location = this->mNodes[index2]->rGetLocation();

    double edge_length = norm_2(GetVectorFromAtoB(node1_location, node2_location));
    return edge_length;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh()
{
    mpDelaunayMesh = NULL;
    this->mMeshChangesDuringSimulation = false;
    Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::~VertexMesh()
{
    Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    /// \todo sort out boundary elements in a vertex mesh (#943)
//    assert(index < this->mBoundaryElements.size() );
    return index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(unsigned elementIndex)
{
    unsigned node_index = UNSIGNED_UNSET;

    if (mVoronoiElementIndexMap.empty())
    {
        node_index = elementIndex;
    }
    else
    {
        for (std::map<unsigned, unsigned>::iterator iter = mVoronoiElementIndexMap.begin();
             iter != mVoronoiElementIndexMap.end();
             ++iter)
        {
            if (iter->second == elementIndex)
            {
                node_index = iter->first;
                break;
            }
        }
    }
    assert(node_index != UNSIGNED_UNSET);
    return node_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(unsigned nodeIndex)
{
    unsigned element_index = UNSIGNED_UNSET;

    if (mVoronoiElementIndexMap.empty())
    {
        element_index = nodeIndex;
    }
    else
    {
        std::map<unsigned, unsigned>::iterator iter = mVoronoiElementIndexMap.find(nodeIndex);

        if (iter == mVoronoiElementIndexMap.end())
        {
            EXCEPTION("This index does not correspond to a VertexElement");
        }
        else
        {
            element_index = iter->second;
        }
    }
    assert(element_index != UNSIGNED_UNSET);
    return element_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // Delete elements
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    mElements.clear();

    // Delete faces
    for (unsigned i=0; i<mFaces.size(); i++)
    {
        delete mFaces[i];
    }
    mFaces.clear();

    // Delete nodes
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
    this->mNodes.clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaces.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM-1, SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetFace(unsigned index) const
{
    assert(index < mFaces.size());
    return mFaces[index];
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCentroidOfElement(unsigned index)
{
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
    unsigned num_nodes_in_element = p_element->GetNumNodes();

    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);

    switch (SPACE_DIM)
    {
        case 1:
        {
            centroid = 0.5*(p_element->GetNodeLocation(0) + p_element->GetNodeLocation(1));
        }
        break;
        case 2: ///\todo Why isn't this just the centre of mass? (#1075)
        {
            c_vector<double, SPACE_DIM> current_node;
            c_vector<double, SPACE_DIM> anticlockwise_node;

            double temp_centroid_x = 0;
            double temp_centroid_y = 0;

            for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
            {
                // Find locations of current node and anticlockwise node
                current_node = p_element->GetNodeLocation(local_index);
                anticlockwise_node = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

                temp_centroid_x += (current_node[0]+anticlockwise_node[0])*(current_node[0]*anticlockwise_node[1]-current_node[1]*anticlockwise_node[0]);
                temp_centroid_y += (current_node[1]+anticlockwise_node[1])*(current_node[0]*anticlockwise_node[1]-current_node[1]*anticlockwise_node[0]);
            }

            double vertex_area = GetVolumeOfElement(index);
            double centroid_coefficient = 1.0/(6.0*vertex_area);

            centroid(0) = centroid_coefficient*temp_centroid_x;
            centroid(1) = centroid_coefficient*temp_centroid_y;
        }
        break;
        case 3:
        {
            for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
            {
                centroid += p_element->GetNodeLocation(local_index);
            }
            centroid /= ((double) num_nodes_in_element);
        }
        break;
        default:
            NEVER_REACHED;
    }
    return centroid;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices;

    // Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = this->GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter = containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
         ++elem_iter)
    {
        // Find the local index of this node in this element
        unsigned local_index = GetElement(*elem_iter)->GetNodeLocalIndex(nodeIndex);

        // Find the global indices of the preceding and successive nodes in this element
        unsigned num_nodes = GetElement(*elem_iter)->GetNumNodes();
        unsigned previous_local_index = (local_index + num_nodes - 1)%num_nodes;
        unsigned next_local_index = (local_index + 1)%num_nodes;

        // Add the global indices of these two nodes to the set of neighbouring node indices
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(previous_local_index));
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(next_local_index));
    }

    return neighbouring_node_indices;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex)
{
    /// \todo We should probably assert here that the node is in fact contained in the element (#1305)

    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices_not_in_this_element;

    // Get the indices of this node's neighbours
    std::set<unsigned> node_neighbours = GetNeighbouringNodeIndices(nodeIndex);

    // Get the indices of the nodes contained in this element
    std::set<unsigned> node_indices_in_this_element;
    for (unsigned local_index=0; local_index<GetElement(elemIndex)->GetNumNodes(); local_index++)
    {
        unsigned global_index = GetElement(elemIndex)->GetNodeGlobalIndex(local_index);
        node_indices_in_this_element.insert(global_index);
    }

    // Check if each neighbour is also in this element; if not, add it to the set
    for (std::set<unsigned>::iterator iter = node_neighbours.begin();
         iter != node_neighbours.end();
         ++iter)
    {
        if (node_indices_in_this_element.find(*iter) == node_indices_in_this_element.end())
        {
            neighbouring_node_indices_not_in_this_element.insert(*iter);
        }
    }

    return neighbouring_node_indices_not_in_this_element;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void VertexMesh<1,1>::ConstructFromMeshReader(AbstractMeshReader<1,1>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void VertexMesh<1,2>::ConstructFromMeshReader(AbstractMeshReader<1,2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void VertexMesh<1,3>::ConstructFromMeshReader(AbstractMeshReader<1,3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void VertexMesh<2,3>::ConstructFromMeshReader(AbstractMeshReader<2,3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void VertexMesh<2,2>::ConstructFromMeshReader(AbstractMeshReader<2,2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned i=0; i<num_nodes; i++)
    {
        node_data = rMeshReader.GetNextNode();
        unsigned is_boundary_node = (unsigned) node_data[2];
        node_data.pop_back();
        this->mNodes.push_back(new Node<2>(i, node_data, is_boundary_node));
    }

    rMeshReader.Reset();

    // Reserve memory for nodes
    mElements.reserve(rMeshReader.GetNumElements());

    // Add elements
    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        // Get the data for this element
        ElementData element_data = rMeshReader.GetNextElementData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j=0; j<num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Use nodes and index to construct this element
        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_index, nodes);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetRegion(attribute_value);
        }
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template<>
void VertexMesh<3,3>::ConstructFromMeshReader(AbstractMeshReader<3,3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned i=0; i<num_nodes; i++)
    {
        node_data = rMeshReader.GetNextNode();
        unsigned is_boundary_node = (unsigned) node_data[3];
        node_data.pop_back();
        this->mNodes.push_back(new Node<3>(i, node_data, is_boundary_node));
    }

    rMeshReader.Reset();

    // Reserve memory for nodes
    mElements.reserve(rMeshReader.GetNumElements());

    // Use a std::set to keep track of which faces have been added to mFaces
    std::set<unsigned> faces_counted;

    // Add elements
    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        ///\todo Horrible hack! (#1076/#1377)
        typedef VertexMeshReader<3,3> VERTEX_MESH_READER;
        assert(dynamic_cast<VERTEX_MESH_READER*>(&rMeshReader) != NULL);

        // Get the data for this element
        VertexElementData element_data = static_cast<VERTEX_MESH_READER*>(&rMeshReader)->GetNextElementDataWithFaces();

        // Get the nodes owned by this element
        std::vector<Node<3>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j=0; j<num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Get the faces owned by this element
        std::vector<VertexElement<2,3>*> faces;
        unsigned num_faces_in_element = element_data.Faces.size();
        for (unsigned i=0; i<num_faces_in_element; i++)
        {
            // Get the data for this face
            ElementData face_data = element_data.Faces[i];

            // Get the face index
            unsigned face_index = face_data.AttributeValue;

            // Get the nodes owned by this face
            std::vector<Node<3>*> nodes_in_face;
            unsigned num_nodes_in_face = face_data.NodeIndices.size();
            for (unsigned j=0; j<num_nodes_in_face; j++)
            {
                assert(face_data.NodeIndices[j] < this->mNodes.size());
                nodes_in_face.push_back(this->mNodes[face_data.NodeIndices[j]]);
            }

            // Use nodes and index to construct this face
            VertexElement<2,3>* p_face = new VertexElement<2,3>(face_index, nodes_in_face);

            // If this face is not already contained in mFaces, add it, and update faces_counted...
            if (faces_counted.find(p_face->GetIndex()) == faces_counted.end())
            {
                mFaces.push_back(p_face);
                faces_counted.insert(p_face->GetIndex());
                faces.push_back(p_face);
            }
            else
            {
                //Don't need the new one.
                delete p_face;
                // ... otherwise use the member of mFaces with this index
                bool face_added = false;
                for (unsigned k=0; k<mFaces.size(); k++)
                {
                    if (mFaces[k]->GetIndex() == face_index)
                    {
                        faces.push_back(mFaces[k]);
                        face_added = true;
                        break;
                    }
                }
                assert(face_added == true);
            }
        }

        ///\todo Store orientations? (#1076/#1377)
        std::vector<bool> orientations = std::vector<bool>(num_faces_in_element, true);

        // Use faces and index to construct this element
        VertexElement<3,3>* p_element = new VertexElement<3,3>(elem_index, faces, orientations, nodes);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetRegion(attribute_value);
        }
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(
    const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector;
    if (mpDelaunayMesh)
    {
        vector = mpDelaunayMesh->GetVectorFromAtoB(rLocationA, rLocationB);
    }
    else
    {
        vector = AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(rLocationA, rLocationB);
    }
    return vector;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVolumeOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3);

    // Get pointer to this element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    double element_volume = 0.0;
    if (SPACE_DIM == 2)
    {
        c_vector<double, SPACE_DIM> first_node = p_element->GetNodeLocation(0);

        unsigned num_nodes_in_element = p_element->GetNumNodes();

        for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
        {
            // Find locations of current node and anticlockwise node
            c_vector<double, SPACE_DIM> current_node = p_element->GetNodeLocation(local_index);
            c_vector<double, SPACE_DIM> anticlockwise_node = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

            /*
             * In order to calculate the area we map the origin to (x[0],y[0])
             * then use GetVectorFromAtoB() to get node cooordiantes
             */
            c_vector<double, SPACE_DIM> transformed_current_node = GetVectorFromAtoB(first_node, current_node);
            c_vector<double, SPACE_DIM> transformed_anticlockwise_node = GetVectorFromAtoB(first_node, anticlockwise_node);

            element_volume += 0.5*(transformed_current_node[0]*transformed_anticlockwise_node[1]
                                   - transformed_anticlockwise_node[0]*transformed_current_node[1]);
        }
    }
    else
    {
        // Loop over faces and add up pyramid volumes
        c_vector<double, SPACE_DIM> pyramid_apex = p_element->GetNodeLocation(0);
        for (unsigned face_index=0; face_index<p_element->GetNumFaces(); face_index++)
        {
            // Get pointer to face
            VertexElement<ELEMENT_DIM-1, SPACE_DIM>* p_face = p_element->GetFace(face_index);

            // Get unit normal to this face
            c_vector<double, SPACE_DIM> unit_normal = GetUnitNormalToFace(p_face);

            // Calculate the perpendicular distance from the plane of the face to the chosen apex
            c_vector<double, SPACE_DIM> base_to_apex = GetVectorFromAtoB(p_face->GetNodeLocation(0), pyramid_apex);
            double perpendicular_distance = inner_prod(base_to_apex, unit_normal);

            // Calculate the area of the face
            double face_area = GetAreaOfFace(p_face);

            // Use these to calculate the volume of the pyramid formed by the face and the point pyramid_apex
            element_volume += face_area * perpendicular_distance / 3;
        }
    }
    return fabs(element_volume);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetSurfaceAreaOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3);

    // Get pointer to this element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    double surface_area = 0.0;
    if (SPACE_DIM == 2)
    {
        unsigned num_nodes_in_element = p_element->GetNumNodes();
        for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
        {
            // Find locations of current node and anticlockwise node
            unsigned current_node_index = p_element->GetNodeGlobalIndex(local_index);
            unsigned anticlockwise_node_index = p_element->GetNodeGlobalIndex((local_index+1)%num_nodes_in_element);

            surface_area += this->GetDistanceBetweenNodes(current_node_index, anticlockwise_node_index);
        }
    }
    else
    {
        // Loop over faces and add up areas
        for (unsigned face_index=0; face_index<p_element->GetNumFaces(); face_index++)
        {
            surface_area += GetAreaOfFace(p_element->GetFace(face_index));
        }
    }
    return surface_area;
}


//////////////////////////////////////////////////////////////////////
//                        2D-specific methods                       //
//////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::ElementIncludesPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);

    // Initialise boolean
    bool element_includes_point = false;

    // Get the element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    // Remap the origin to the first vertex to allow alternative distance metrics to be used in subclasses
    c_vector<double, SPACE_DIM> first_vertex = p_element->GetNodeLocation(0);

    c_vector<double, SPACE_DIM> test_point =  GetVectorFromAtoB(first_vertex,rTestPoint);

    // Loop over edges of the element
    for (unsigned local_index=0; local_index<num_nodes; local_index++)
    {
        // Get the end points of this edge
        // Remap to the origin to allow alternative distance metrics to be used in subclasses
        c_vector<double, SPACE_DIM> vertexA = GetVectorFromAtoB(first_vertex, p_element->GetNodeLocation(local_index));
        c_vector<double, SPACE_DIM> vertexB = GetVectorFromAtoB(first_vertex, p_element->GetNodeLocation((local_index+1)%num_nodes));


        // Check if this edge crosses the ray running out horizontally (increasing x, fixed y) from the test point
        c_vector<double, SPACE_DIM> vector_a_to_point = GetVectorFromAtoB(vertexA, test_point);
        c_vector<double, SPACE_DIM> vector_b_to_point = GetVectorFromAtoB(vertexB, test_point);
        c_vector<double, SPACE_DIM> vector_a_to_b = GetVectorFromAtoB(vertexA, vertexB);

        // Pathological case - test point coincides with vertexA or vertexB
        if (    (norm_2(vector_a_to_point) < DBL_EPSILON)
             || (norm_2(vector_b_to_point) < DBL_EPSILON) )
        {
            return false;
        }

        // Pathological case - ray coincides with horizontal edge
        if ( (fabs(vector_a_to_b[1]) < DBL_EPSILON) &&
             (fabs(vector_a_to_point[1]) < DBL_EPSILON) &&
             (fabs(vector_b_to_point[1]) < DBL_EPSILON) )
        {
            if ( (vector_a_to_point[0]>0) != (vector_b_to_point[0]>0) )
            {
                return false;
            }
        }

        // Non-pathological case
        // A and B on different sides of the line y = test_point[1]
        if ( (vertexA[1] > test_point[1]) != (vertexB[1] > test_point[1]) )
        {
            // intersection of y=test_point[1] and vector_a_to_b is on the Right of test_point
            if (test_point[0] < vertexA[0] + vector_a_to_b[0]*vector_a_to_point[1]/vector_a_to_b[1])
            {
                element_includes_point = !element_includes_point;
            }
        }
    }
    return element_includes_point;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetLocalIndexForElementEdgeClosestToPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);

    // Get the element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    double min_squared_normal_distance = DBL_MAX;
    unsigned min_distance_edge_index = UINT_MAX;

    // Loop over edges of the element
    for (unsigned local_index=0; local_index<num_nodes; local_index++)
    {
        // Get the end points of this edge
        c_vector<double, SPACE_DIM> vertexA = p_element->GetNodeLocation(local_index);
        c_vector<double, SPACE_DIM> vertexB = p_element->GetNodeLocation((local_index+1)%num_nodes);

        c_vector<double, SPACE_DIM> vector_a_to_point = this->GetVectorFromAtoB(vertexA, rTestPoint);
        c_vector<double, SPACE_DIM> vector_a_to_b = this->GetVectorFromAtoB(vertexA, vertexB);

        c_vector<double, SPACE_DIM> edge_ab_unit_vector = vector_a_to_b/norm_2(vector_a_to_b);
        double distance_parallel_to_edge = inner_prod(vector_a_to_point, edge_ab_unit_vector);

        double squared_distance_normal_to_edge = SmallPow(norm_2(vector_a_to_point), 2) - SmallPow(distance_parallel_to_edge, 2);

        // Make sure node is within the confines of the edge and is the nearest edge to the node \this breks for convex elements
        if (squared_distance_normal_to_edge < min_squared_normal_distance &&
            distance_parallel_to_edge >=0 &&
            distance_parallel_to_edge <= norm_2(vector_a_to_b))
        {
            min_squared_normal_distance = squared_distance_normal_to_edge;
            min_distance_edge_index = local_index;
        }
    }
    return min_distance_edge_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> VertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMomentsOfElement(unsigned index)
{
    assert(SPACE_DIM == 2);

    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
    unsigned num_nodes_in_element = p_element->GetNumNodes();
    c_vector<double, 2> centroid = GetCentroidOfElement(index);

    c_vector<double, 3> moments = zero_vector<double>(3);

    unsigned node_1;
    unsigned node_2;

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        node_1 = local_index;
        node_2 = (local_index+1)%num_nodes_in_element;

        // Original position of nodes
        c_vector<double, 2> original_pos_1 = p_element->GetNodeLocation(node_1);
        c_vector<double, 2> original_pos_2 = p_element->GetNodeLocation(node_2);

        // Node position so centerd on origin
        c_vector<double, 2> pos_1 = this->GetVectorFromAtoB(centroid, original_pos_1);
        c_vector<double, 2> pos_2 = this->GetVectorFromAtoB(centroid, original_pos_2);

        /*
         * Note these formulae require the polygon to be centered on the origin
         */
        double a = pos_1(0)*pos_2(1)-pos_2(0)*pos_1(1);

        // Ixx
        moments(0) += (  pos_1(1)*pos_1(1)
                       + pos_1(1)*pos_2(1)
                       + pos_2(1)*pos_2(1) ) * a;

        // Iyy
        moments(1) += (  pos_1(0)*pos_1(0)
                       + pos_1(0)*pos_2(0)
                       + pos_2(0)*pos_2(0) ) * a;

        // Ixy
        moments(2) += (  pos_1(0)*pos_2(1)
                       + 2*pos_1(0)*pos_1(1)
                       + 2*pos_2(0)*pos_2(1)
                       + pos_2(0)*pos_1(1) ) * a;
    }

    moments(0) /= 12;
    moments(1) /= 12;
    moments(2) /= 24;

    return moments;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetShortAxisOfElement(unsigned index)
{
    assert(SPACE_DIM == 2);

    c_vector<double, SPACE_DIM> short_axis = zero_vector<double>(SPACE_DIM);
    c_vector<double, 3> moments = CalculateMomentsOfElement(index);

    double largest_eigenvalue, discriminant;

    discriminant = sqrt((moments(0) - moments(1))*(moments(0) - moments(1)) + 4.0*moments(2)*moments(2));
    // This is always the largest eigenvalue as both eigenvalues are real as it is a
    // symmetric matrix
    largest_eigenvalue = (moments(0) + moments(1) + discriminant)*0.5;
    if (fabs(discriminant) < 1e-10)
    {
        // Return a random unit vector
        short_axis(0) = RandomNumberGenerator::Instance()->ranf();
        short_axis(1) = sqrt(1.0 - short_axis(0)*short_axis(0));
    }
    else
    {
        if (moments(2) == 0.0)
        {
            short_axis(0) = (moments(0) < moments(1)) ? 0.0 : 1.0;
            short_axis(1) = (moments(0) < moments(1)) ? 1.0 : 0.0;
        }
        else
        {
            short_axis(0) = 1.0;
            short_axis(1) = (moments(0) - largest_eigenvalue)/moments(2);

            double length_short_axis = norm_2(short_axis);

            short_axis /= length_short_axis;
        }
    }
    return short_axis;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetAreaGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex)
{
    assert(SPACE_DIM==2);

    unsigned num_nodes_in_element = pElement->GetNumNodes();
    unsigned next_local_index = (localIndex+1)%num_nodes_in_element;

    // We add an extra num_nodes_in_element in the line below as otherwise this term can be negative, which breaks the % operator
    unsigned previous_local_index = (num_nodes_in_element+localIndex-1)%num_nodes_in_element;

    c_vector<double, SPACE_DIM> previous_node_location = pElement->GetNodeLocation(previous_local_index);
    c_vector<double, SPACE_DIM> next_node_location = pElement->GetNodeLocation(next_local_index);
    c_vector<double, SPACE_DIM> difference_vector = this->GetVectorFromAtoB(previous_node_location, next_node_location);

    c_vector<double, SPACE_DIM> area_gradient;

    area_gradient[0] = 0.5*difference_vector[1];
    area_gradient[1] = -0.5*difference_vector[0];

    return area_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetPreviousEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex)
{
    assert(SPACE_DIM==2);

    unsigned num_nodes_in_element = pElement->GetNumNodes();

    // We add an extra localIndex-1 in the line below as otherwise this term can be negative, which breaks the % operator
    unsigned previous_local_index = (num_nodes_in_element+localIndex-1)%num_nodes_in_element;

    unsigned current_global_index = pElement->GetNodeGlobalIndex(localIndex);
    unsigned previous_global_index = pElement->GetNodeGlobalIndex(previous_local_index);

    double previous_edge_length = this->GetDistanceBetweenNodes(current_global_index, previous_global_index);
    assert(previous_edge_length > DBL_EPSILON);

    c_vector<double, SPACE_DIM> previous_edge_gradient = this->GetVectorFromAtoB(pElement->GetNodeLocation(previous_local_index), pElement->GetNodeLocation(localIndex))/previous_edge_length;

    return previous_edge_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNextEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex)
{
    assert(SPACE_DIM==2);

    unsigned next_local_index = (localIndex+1)%(pElement->GetNumNodes());

    unsigned current_global_index = pElement->GetNodeGlobalIndex(localIndex);
    unsigned next_global_index = pElement->GetNodeGlobalIndex(next_local_index);

    double next_edge_length = this->GetDistanceBetweenNodes(current_global_index, next_global_index);
    assert(next_edge_length > DBL_EPSILON);

    c_vector<double, SPACE_DIM> next_edge_gradient = this->GetVectorFromAtoB(pElement->GetNodeLocation(next_local_index), pElement->GetNodeLocation(localIndex))/next_edge_length;

    return next_edge_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetPerimeterGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex)
{
    assert(SPACE_DIM==2);

    c_vector<double, SPACE_DIM> previous_edge_gradient = GetPreviousEdgeGradientOfElementAtNode(pElement, localIndex);
    c_vector<double, SPACE_DIM> next_edge_gradient = GetNextEdgeGradientOfElementAtNode(pElement, localIndex);

    return previous_edge_gradient + next_edge_gradient;
}


//////////////////////////////////////////////////////////////////////
//                        3D-specific methods                       //
//////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetUnitNormalToFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace)
{
    assert(SPACE_DIM == 3);

    // As we are in 3D, the face must have at least three vertices, so use its first three vertices
    c_vector<double, SPACE_DIM> v0 = pFace->GetNode(0)->rGetLocation();
    c_vector<double, SPACE_DIM> v1 = pFace->GetNode(1)->rGetLocation();
    c_vector<double, SPACE_DIM> v2 = pFace->GetNode(2)->rGetLocation();

    c_vector<double, SPACE_DIM> v1_minus_v0 = this->GetVectorFromAtoB(v0, v1);
    c_vector<double, SPACE_DIM> v2_minus_v0 = this->GetVectorFromAtoB(v0, v2);

    c_vector<double, SPACE_DIM> unit_normal = zero_vector<double>(SPACE_DIM);
    unit_normal(0) = v1_minus_v0(1)*v2_minus_v0(2) - v1_minus_v0(2)*v2_minus_v0(1);
    unit_normal(1) = v1_minus_v0(2)*v2_minus_v0(0) - v1_minus_v0(0)*v2_minus_v0(2);
    unit_normal(2) = v1_minus_v0(0)*v2_minus_v0(1) - v1_minus_v0(1)*v2_minus_v0(0);

    // Normalize the normal vector
    unit_normal /= norm_2(unit_normal);

    return unit_normal;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetAreaOfFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace)
{
    assert(SPACE_DIM == 3);

    // Get the unit normal to the plane of this face
    c_vector<double, SPACE_DIM> unit_normal = GetUnitNormalToFace(pFace);

    // Select the largest absolute coordinate to ignore for planar projection
    double abs_x = unit_normal[0]>0 ? unit_normal[0]>0 : -unit_normal[0];
    double abs_y = unit_normal[1]>0 ? unit_normal[1]>0 : -unit_normal[1];
    double abs_z = unit_normal[2]>0 ? unit_normal[2]>0 : -unit_normal[2];

    unsigned dim_to_ignore = 2; // ignore z coordinate
    double abs = abs_z;

    if (abs_x > abs_y)
    {
        if (abs_x > abs_z)
        {
            dim_to_ignore = 0; // ignore x coordinate
            abs = abs_x;
        }
    }
    else if (abs_y > abs_z)
    {
        dim_to_ignore = 1; // ignore y coordinate
        abs = abs_y;
    }

    // Compute area of the 2D projection
    c_vector<double, SPACE_DIM-1> current_vertex;
    c_vector<double, SPACE_DIM-1> anticlockwise_vertex;

    unsigned num_nodes_in_face = pFace->GetNumNodes();

    unsigned dim1 = dim_to_ignore==0 ? 1 : 0;
    unsigned dim2 = dim_to_ignore==2 ? 1 : 2;

    double face_area = 0.0;
    for (unsigned local_index=0; local_index<num_nodes_in_face; local_index++)
    {
        // Find locations of current vertex and anticlockwise vertex
        current_vertex[0] = pFace->GetNodeLocation(local_index, dim1);
        current_vertex[1] = pFace->GetNodeLocation(local_index, dim2);
        anticlockwise_vertex[0] = pFace->GetNodeLocation((local_index+1)%num_nodes_in_face, dim1);
        anticlockwise_vertex[1] = pFace->GetNodeLocation((local_index+1)%num_nodes_in_face, dim2);

        // It doesn't matter if the face is oriented clockwise or not, since we area only interested in the area
        face_area += 0.5*(current_vertex[0]*anticlockwise_vertex[1] - anticlockwise_vertex[0]*current_vertex[1]);
    }

    // Scale to get area before projection
    face_area /= abs;
    return fabs(face_area);
}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class VertexMesh<1,1>;
template class VertexMesh<1,2>;
template class VertexMesh<1,3>;
template class VertexMesh<2,2>;
template class VertexMesh<2,3>;
template class VertexMesh<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMesh)
