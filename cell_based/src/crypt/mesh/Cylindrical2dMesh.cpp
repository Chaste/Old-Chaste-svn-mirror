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

#include <map>

/* These lines are very useful for debugging (visualize with 'showme').
#include "TrianglesMeshWriter.hpp"
TrianglesMeshWriter<2,2> mesh_writer("Cylindrical2dMeshDebug", "mesh", false);
mesh_writer.WriteFilesUsingMesh(*this);
*/
#include "Cylindrical2dMesh.hpp"


Cylindrical2dMesh::Cylindrical2dMesh(double width)
  : MutableMesh<2,2>(),
    mWidth(width)
{
    assert(width > 0.0);
}


Cylindrical2dMesh::~Cylindrical2dMesh()
{
}


Cylindrical2dMesh::Cylindrical2dMesh(double width, std::vector<Node<2>* > nodes)
  : MutableMesh<2,2>(),
    mWidth(width)
{
    assert(width > 0.0);
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<2>* p_temp_node = nodes[index];
        double x = p_temp_node->rGetLocation()[0];
        x = x; // Fix optimised build
        assert( 0 <= x && x < width);
        mNodes.push_back(p_temp_node);
    }

    NodeMap node_map(nodes.size());
    ReMesh(node_map);
}


void Cylindrical2dMesh::UpdateTopAndBottom()
{
    ChasteCuboid<2> bounding_box=CalculateBoundingBox();
    mBottom = bounding_box.rGetLowerCorner()[1];
    mTop = bounding_box.rGetUpperCorner()[1];
}


void Cylindrical2dMesh::CreateMirrorNodes()
{
    double half_way = 0.5*mWidth;

    mLeftOriginals.clear();
    mLeftImages.clear();
    mImageToLeftOriginalNodeMap.clear();
    mRightOriginals.clear();
    mRightImages.clear();
    mImageToRightOriginalNodeMap.clear();
    mLeftPeriodicBoundaryElementIndices.clear();
    mRightPeriodicBoundaryElementIndices.clear();

    for (AbstractMesh<2,2>::NodeIterator node_iter = GetNodeIteratorBegin();
         node_iter != GetNodeIteratorEnd();
         ++node_iter)
    {
        c_vector<double, 2> location = node_iter->rGetLocation();
        unsigned this_node_index = node_iter->GetIndex();
        double this_node_x_location = location[0];

        // Check the mesh currently conforms to the dimensions given
        assert(0.0 <= location[0]);
        assert(location[0] <= mWidth);

        // Put the nodes which are to be mirrored in the relevant vectors
        if (this_node_x_location < half_way)
        {
            mLeftOriginals.push_back(this_node_index);
        }
        else
        {
            mRightOriginals.push_back(this_node_index);
        }
    }

    // Go through the left original nodes and create an image node
    // recording its new index
    for (unsigned i=0; i<mLeftOriginals.size(); i++)
    {
        c_vector<double, 2> location = mNodes[mLeftOriginals[i]]->rGetLocation();
        location[0] = location[0] + mWidth;

        unsigned new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));
        mLeftImages.push_back(new_node_index);
        mImageToLeftOriginalNodeMap[new_node_index] = mLeftOriginals[i];
    }

    // Go through the right original nodes and create an image node
    // recording its new index
    for (unsigned i=0; i<mRightOriginals.size(); i++)
    {
        // Create new image nodes
        c_vector<double, 2> location = mNodes[mRightOriginals[i]]->rGetLocation();
        location[0] = location[0] - mWidth;

        unsigned new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));
        mRightImages.push_back(new_node_index);
        mImageToRightOriginalNodeMap[new_node_index] = mRightOriginals[i];
    }

    assert(mRightOriginals.size()==mRightImages.size());
    assert(mLeftOriginals.size()==mLeftImages.size());
    assert(mImageToLeftOriginalNodeMap.size()==mLeftOriginals.size());
    assert(mImageToRightOriginalNodeMap.size()==mRightOriginals.size());
}


void Cylindrical2dMesh::CreateHaloNodes()
{
    UpdateTopAndBottom();

    mTopHaloNodes.clear();
    mBottomHaloNodes.clear();

    unsigned num_halo_nodes = (unsigned)(floor(mWidth*2.0));
    double halo_node_separation = mWidth/((double)(num_halo_nodes));
    double y_top_coordinate = mTop + halo_node_separation;
    double y_bottom_coordinate = mBottom - halo_node_separation;

    c_vector<double, 2> location;
    for (unsigned i=0; i< num_halo_nodes; i++)
    {
       double x_coordinate = 0.5*halo_node_separation + (double)(i)*halo_node_separation;

       // Inserting top halo node in mesh
       location[0] = x_coordinate;
       location[1] = y_top_coordinate;
       unsigned new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));
       mTopHaloNodes.push_back(new_node_index);

       location[1] = y_bottom_coordinate;
       new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));
       mBottomHaloNodes.push_back(new_node_index);
    }
}


void Cylindrical2dMesh::ReMesh(NodeMap &map)
{
    unsigned old_num_all_nodes = GetNumAllNodes();

    map.Resize(old_num_all_nodes);
    map.ResetToIdentity();

    // Flag the deleted nodes as deleted in the map
    for (unsigned i=0; i<old_num_all_nodes; i++)
    {
        if (mNodes[i]->IsDeleted())
        {
            map.SetDeleted(i);
        }
    }

    CreateHaloNodes();

    // Create a mirrored load of nodes for the normal remesher to work with.
    CreateMirrorNodes();

    /*
     * The mesh now has messed up boundary elements, but this
     * doesn't matter as the ReMesh below doesn't read them in
     * and reconstructs the boundary elements.
     */

    // Call the normal re-mesh. Note that the mesh now has lots
    // of extra nodes which will be deleted, hence the name 'big_map'
    NodeMap big_map(GetNumAllNodes());
    MutableMesh<2,2>::ReMesh(big_map);

    /*
     * If the big_map isn't the identity map, the little map ('map') needs to be
     * altered accordingly before being passed to the user. Not sure how this all
     * works, so deal with this bridge when we get to it.
     */
    assert(big_map.IsIdentityMap());

    // Re-index the vectors according to the big nodemap, and set up the maps.
    mImageToLeftOriginalNodeMap.clear();
    mImageToRightOriginalNodeMap.clear();

    assert(mLeftOriginals.size()==mLeftImages.size());
    assert(mRightOriginals.size()==mRightImages.size());

    for (unsigned i=0; i<mLeftOriginals.size(); i++)
    {
        mLeftOriginals[i] = big_map.GetNewIndex(mLeftOriginals[i]);
        mLeftImages[i] = big_map.GetNewIndex(mLeftImages[i]);
        mImageToLeftOriginalNodeMap[mLeftImages[i]] = mLeftOriginals[i];
    }

    for (unsigned i=0; i<mRightOriginals.size(); i++)
    {
        mRightOriginals[i] = big_map.GetNewIndex(mRightOriginals[i]);
        mRightImages[i] = big_map.GetNewIndex(mRightImages[i]);
        mImageToRightOriginalNodeMap[mRightImages[i]] = mRightOriginals[i];
    }

    for (unsigned i=0; i<mTopHaloNodes.size(); i++)
    {
        mTopHaloNodes[i] = big_map.GetNewIndex(mTopHaloNodes[i]);
        mBottomHaloNodes[i] = big_map.GetNewIndex(mBottomHaloNodes[i]);
    }


    // This method checks elements crossing the periodic boundary have been meshed in the same way at each side.
    CorrectNonPeriodicMesh();

    /*
     * This method takes in the double-sized mesh, with its new boundary
     * elements, and removes the relevant nodes, elements and boundary
     * elements to leave a proper periodic mesh.
     */
    ReconstructCylindricalMesh();

    DeleteHaloNodes();

    /*
     * Create a random boundary element between two nodes of the first
     * element if it is not deleted. This is a temporary measure to get
     * around re-index crashing when there are no boundary elements.
     */
    unsigned num_elements = GetNumAllElements();
    bool boundary_element_made = false;
    unsigned elem_index = 0;

    while (elem_index<num_elements && !boundary_element_made)
    {
        Element<2,2>* p_element = GetElement(elem_index);
        if (!p_element->IsDeleted())
        {
            boundary_element_made = true;
            std::vector<Node<2>*> nodes;
            nodes.push_back(p_element->GetNode(0));
            nodes.push_back(p_element->GetNode(1));
            BoundaryElement<1,2>* p_boundary_element = new BoundaryElement<1,2>(0, nodes);
            p_boundary_element->RegisterWithNodes();
            mBoundaryElements.push_back(p_boundary_element);
            this->mBoundaryElementWeightedDirections.push_back(zero_vector<double>(2));
            this->mBoundaryElementJacobianDeterminants.push_back(0.0);
        }
        elem_index++;
    }

    // Now call ReIndex to remove the temporary nodes which are marked as deleted.
    NodeMap reindex_map(GetNumAllNodes());
    ReIndex(reindex_map);
    assert(!reindex_map.IsIdentityMap());  // maybe don't need this

    // Go through the reindex map and use it to populate the original NodeMap
    // (the one that is returned to the user)
    for (unsigned i=0; i<map.Size(); i++) // only going up to be size of map, not size of reindex_map
    {
        if (reindex_map.IsDeleted(i))
        {
            /*
             * i < num_original_nodes and node is deleted, this should correspond to
             * a node that was labelled as before the remeshing, so should have already
             * been set as deleted in the map above.
             */
            assert(map.IsDeleted(i));
        }
        else
        {
            map.SetNewIndex(i, reindex_map.GetNewIndex(i) );
        }
    }

    // We can now clear the index vectors & maps; they are only used for remeshing
    mLeftOriginals.clear();
    mLeftImages.clear();
    mImageToLeftOriginalNodeMap.clear();
    mRightOriginals.clear();
    mRightImages.clear();
    mImageToRightOriginalNodeMap.clear();
    mLeftPeriodicBoundaryElementIndices.clear();
    mRightPeriodicBoundaryElementIndices.clear();
}


void Cylindrical2dMesh::ReconstructCylindricalMesh()
{
    /*
     * Figure out which elements have real nodes and image nodes in them
     * and replace image nodes with corresponding real ones.
     */
    for (MutableMesh<2,2>::ElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        // Left images are on the right of the mesh
        unsigned number_of_left_image_nodes = 0;
        unsigned number_of_right_image_nodes = 0;

        for (unsigned i=0; i<3; i++)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);

            if (mImageToLeftOriginalNodeMap.find(this_node_index) != mImageToLeftOriginalNodeMap.end())
            {
                number_of_left_image_nodes++;
            }
            else if (mImageToRightOriginalNodeMap.find(this_node_index) != mImageToRightOriginalNodeMap.end())
            {
                number_of_right_image_nodes++;
            }
        }

        // Delete all the elements on the left hand side (images of right)...
        if (number_of_right_image_nodes >= 1)
        {
            elem_iter->MarkAsDeleted();
            mDeletedElementIndices.push_back(elem_iter->GetIndex());
        }

        // Delete only purely imaginary elements on the right (images of left nodes)
        if (number_of_left_image_nodes == 3)
        {
            elem_iter->MarkAsDeleted();
            mDeletedElementIndices.push_back(elem_iter->GetIndex());
        }

        /*
         * If some are images then replace them with the real nodes. There
         * can be elements with either two image nodes on the right (and one
         * real) or one image node on the right (and two real).
         */
        if (number_of_left_image_nodes == 1 || number_of_left_image_nodes == 2)
        {
            for (unsigned i=0; i<3; i++)
            {
                unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);
                std::map<unsigned, unsigned>::iterator it = mImageToLeftOriginalNodeMap.find(this_node_index);
                if (it != mImageToLeftOriginalNodeMap.end())
                {
                    elem_iter->ReplaceNode(mNodes[this_node_index], mNodes[it->second]);
                }
            }
        }
    } // end of loop over elements

    /*
     * Figure out which boundary elements have real nodes and image nodes in them
     * and replace image nodes with corresponding real ones.
     */
    for (unsigned elem_index = 0; elem_index<GetNumAllBoundaryElements(); elem_index++)
    {
        BoundaryElement<1,2>* p_boundary_element = GetBoundaryElement(elem_index);
        if (!p_boundary_element->IsDeleted())
        {
            unsigned number_of_image_nodes = 0;
            for (unsigned i=0; i<2; i++)
            {
                unsigned this_node_index = p_boundary_element->GetNodeGlobalIndex(i);

                if (mImageToLeftOriginalNodeMap.find(this_node_index)
                    != mImageToLeftOriginalNodeMap.end())
                {
                    number_of_image_nodes++;
                }
                else if (mImageToRightOriginalNodeMap.find(this_node_index)
                    != mImageToRightOriginalNodeMap.end())
                {
                    number_of_image_nodes++;
                }
            }

            if (number_of_image_nodes==2)
            {
                p_boundary_element->MarkAsDeleted();
                mDeletedBoundaryElementIndices.push_back(p_boundary_element->GetIndex());
            }

            /*
             * To avoid having two copies of the boundary elements on the periodic
             * boundaries we only deal with the elements on the left image and
             * delete the ones on the right image.
             */
            if (number_of_image_nodes==1)
            {
                for (unsigned i=0; i<2; i++)
                {
                    unsigned this_node_index = p_boundary_element->GetNodeGlobalIndex(i);
                    std::map<unsigned, unsigned>::iterator it = mImageToLeftOriginalNodeMap.find(this_node_index);
                    if (it != mImageToLeftOriginalNodeMap.end())
                    {
                        p_boundary_element->ReplaceNode(mNodes[this_node_index], mNodes[it->second]);
                    }
                    else
                    {
                        it = mImageToRightOriginalNodeMap.find(this_node_index);
                        if (it != mImageToRightOriginalNodeMap.end())
                        {
                            p_boundary_element->MarkAsDeleted();
                            mDeletedBoundaryElementIndices.push_back(p_boundary_element->GetIndex());
                        }
                    }
                }
            }
        }
    }

    // Delete all image nodes unless they have already gone (halo nodes)
    for (unsigned i=0; i<mLeftImages.size(); i++)
    {
        mNodes[mLeftImages[i]]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(mLeftImages[i]);
    }

    for (unsigned i=0; i<mRightImages.size(); i++)
    {
        mNodes[mRightImages[i]]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(mRightImages[i]);
    }
}


void Cylindrical2dMesh::DeleteHaloNodes()
{
    assert(mTopHaloNodes.size()==mBottomHaloNodes.size());
    for (unsigned i=0; i<mTopHaloNodes.size(); i++)
    {
        DeleteBoundaryNodeAt(mTopHaloNodes[i]);
        DeleteBoundaryNodeAt(mBottomHaloNodes[i]);
    }
}


c_vector<double, 2> Cylindrical2dMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth > 0.0);

    c_vector<double, 2> location1 = rLocation1;
    c_vector<double, 2> location2 = rLocation2;

    location1[0] = fmod(location1[0], mWidth);
    location2[0] = fmod(location2[0], mWidth);

    c_vector<double, 2> vector = location2 - location1;

    /*
     * Handle the cylindrical condition here: if the points are more
     * than halfway around the cylinder apart, measure the other way.
     */
    if (vector[0] > 0.5*mWidth)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -0.5*mWidth)
    {
        vector[0] += mWidth;
    }
    return vector;
}


void Cylindrical2dMesh::SetNode(unsigned index, ChastePoint<2> point, bool concreteMove)
{
    // Perform a periodic movement if necessary
    if (point.rGetLocation()[0] >= mWidth)
    {
        // Move point to the left
        point.SetCoordinate(0, point.rGetLocation()[0] - mWidth);
    }
    else if (point.rGetLocation()[0] < 0.0)
    {
        // Move point to the right
        point.SetCoordinate(0, point.rGetLocation()[0] + mWidth);
    }

    // Update the node's location
    MutableMesh<2,2>::SetNode(index, point, concreteMove);
}


double Cylindrical2dMesh::GetWidth(const unsigned& rDimension) const
{
    double width = 0.0;
    assert(rDimension==0 || rDimension==1);
    if (rDimension==0)
    {
        width = mWidth;
    }
    else
    {
        width = MutableMesh<2,2>::GetWidth(rDimension);
    }
    return width;
}


unsigned Cylindrical2dMesh::AddNode(Node<2>* pNewNode)
{
    unsigned node_index = MutableMesh<2,2>::AddNode(pNewNode);

    // If necessary move it to be back on the cylinder
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point, false);

    return node_index;
}


void Cylindrical2dMesh::CorrectNonPeriodicMesh()
{
    GenerateVectorsOfElementsStraddlingPeriodicBoundaries();

    /*
     * Copy the member variables into new vectors, which we modify
     * by knocking out elements which pair up on each side.
     */
    std::set<unsigned> temp_left_hand_side_elements = mLeftPeriodicBoundaryElementIndices;
    std::set<unsigned> temp_right_hand_side_elements = mRightPeriodicBoundaryElementIndices;
    assert(mLeftPeriodicBoundaryElementIndices.size()==mRightPeriodicBoundaryElementIndices.size());

    // Go through all of the elements on the left periodic boundary
    for (std::set<unsigned>::iterator left_iter = mLeftPeriodicBoundaryElementIndices.begin();
         left_iter != mLeftPeriodicBoundaryElementIndices.end();
         ++left_iter)
    {
        unsigned elem_index = *left_iter;
        Element<2,2>* p_element = GetElement(elem_index);

        /*
         * Make lists of the nodes which the elements on the left contain and
         * the nodes which should be in a corresponding element on the right.
         */
        c_vector<unsigned,3> original_element_node_indices;
        c_vector<unsigned,3> corresponding_element_node_indices;
        for (unsigned i=0; i<3; i++)
        {
            original_element_node_indices[i] = p_element->GetNodeGlobalIndex(i);
            corresponding_element_node_indices[i] = GetCorrespondingNodeIndex(original_element_node_indices[i]);
        }

        // Search the right hand side elements for the corresponding element
        for (std::set<unsigned>::iterator right_iter = mRightPeriodicBoundaryElementIndices.begin();
             right_iter != mRightPeriodicBoundaryElementIndices.end();
             ++right_iter)
        {
            unsigned corresponding_elem_index = *right_iter;
            Element<2,2>* p_corresponding_element = GetElement(corresponding_elem_index);

            bool is_corresponding_node = true;

            for (unsigned i=0; i<3; i++)
            {
                if ( (corresponding_element_node_indices[i] != p_corresponding_element->GetNodeGlobalIndex(0)) &&
                     (corresponding_element_node_indices[i] != p_corresponding_element->GetNodeGlobalIndex(1)) &&
                     (corresponding_element_node_indices[i] != p_corresponding_element->GetNodeGlobalIndex(2)) )
                {
                    is_corresponding_node = false;
                    break;
                }
            }

            if (is_corresponding_node)
            {
                // Remove original and corresponding element from sets
                temp_left_hand_side_elements.erase(elem_index);
                temp_right_hand_side_elements.erase(corresponding_elem_index);
            }
        }
    }

    /*
     * If either of these ever throw you have more than one situation where the mesher has an option
     * of how to mesh. If it does ever throw you need to be cleverer and match up the
     * elements into as many pairs as possible on the left hand and right hand sides.
     */
    assert(temp_left_hand_side_elements.size() <= 2);
    assert(temp_right_hand_side_elements.size() <= 2);

    /*
     * Now we just have to use the first pair of elements and copy their info over to the other side.
     * First we need to get hold of both elements on either side.
     */
    if (temp_left_hand_side_elements.empty() || temp_right_hand_side_elements.empty())
    {
        assert(temp_right_hand_side_elements.empty());
        assert(temp_left_hand_side_elements.empty());
    }
    else
    {
        assert(temp_right_hand_side_elements.size() == 2 && temp_left_hand_side_elements.size() == 2);
        if (temp_right_hand_side_elements.size() == 2)
        {
            // Use the right hand side meshing and map to left
            UseTheseElementsToDecideMeshing(temp_right_hand_side_elements);
        }
        else
        {
            /*
             * If you get here there are more than two mixed up elements on the periodic edge.
             * We need to knock the pair out and then rerun this function. This shouldn't be
             * too hard to do but is as yet unnecessary.
             */
            NEVER_REACHED;
        }
    }
}

void Cylindrical2dMesh::UseTheseElementsToDecideMeshing(std::set<unsigned>& rMainSideElements)
{
    assert(rMainSideElements.size() == 2);

    // We find the four nodes surrounding the dodgy meshing, on each side.
    std::set<unsigned> main_four_nodes;
    for (std::set<unsigned>::iterator left_iter = rMainSideElements.begin();
         left_iter != rMainSideElements.end();
         ++left_iter)
    {
        unsigned elem_index = *left_iter;
        Element<2,2>* p_element = GetElement(elem_index);
        for (unsigned i=0; i<3; i++)
        {
            unsigned index = p_element->GetNodeGlobalIndex(i);
            main_four_nodes.insert(index);
        }
    }
    assert(main_four_nodes.size() == 4);

    std::set<unsigned> other_four_nodes;
    for (std::set<unsigned>::iterator iter = main_four_nodes.begin();
         iter != main_four_nodes.end();
         ++iter)
    {
        other_four_nodes.insert(GetCorrespondingNodeIndex(*iter));
    }
    assert(other_four_nodes.size() == 4);

    /*
     * Find the elements surrounded by the nodes on the right
     * and change them to match the elements on the left.
     */
    std::vector<unsigned> corresponding_elements;

    // Loop over all elements
    for (MutableMesh<2,2>::ElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        // Loop over the nodes of the element
        if (!(other_four_nodes.find(elem_iter->GetNodeGlobalIndex(0))==other_four_nodes.end()) &&
            !(other_four_nodes.find(elem_iter->GetNodeGlobalIndex(1))==other_four_nodes.end()) &&
            !(other_four_nodes.find(elem_iter->GetNodeGlobalIndex(2))==other_four_nodes.end()) )
        {
            corresponding_elements.push_back(elem_iter->GetIndex());
            elem_iter->MarkAsDeleted();
            mDeletedElementIndices.push_back(elem_iter->GetIndex());
        }
    }
    assert(corresponding_elements.size() == 2);

    // Now corresponding_elements contains the two elements which are going to be replaced by rMainSideElements
    unsigned num_elements = GetNumAllElements();
    for (std::set<unsigned>::iterator iter = rMainSideElements.begin();
         iter != rMainSideElements.end();
         ++iter)
    {
        Element<2,2>* p_main_element = GetElement(*iter);
        std::vector<Node<2>*> nodes;

        // Put corresponding nodes into a std::vector
        for (unsigned i=0; i<3; i++)
        {
            unsigned main_node = p_main_element->GetNodeGlobalIndex(i);
            nodes.push_back(this->GetNode(GetCorrespondingNodeIndex(main_node)));
        }

        // Make a new element
        Element<2,2>* p_new_element = new Element<2,2>(num_elements, nodes);
        this->mElements.push_back(p_new_element);
        this->mElementJacobians.push_back(zero_matrix<double>(2,2));
        this->mElementInverseJacobians.push_back(zero_matrix<double>(2,2));
        this->mElementJacobianDeterminants.push_back(0.0);
        num_elements++;
    }

    // Reindex to get rid of extra elements indices
    NodeMap map(GetNumAllNodes());
    this->ReIndex(map);
}


void Cylindrical2dMesh::GenerateVectorsOfElementsStraddlingPeriodicBoundaries()
{
    mLeftPeriodicBoundaryElementIndices.clear();
    mRightPeriodicBoundaryElementIndices.clear();

    for (MutableMesh<2,2>::ElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        // Left images are on the right of the mesh
        unsigned number_of_left_image_nodes = 0;
        unsigned number_of_right_image_nodes = 0;
        for (unsigned i=0; i<3; i++)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);
            if (mImageToLeftOriginalNodeMap.find(this_node_index) != mImageToLeftOriginalNodeMap.end())
            {
                number_of_left_image_nodes++;
            }
            else if (mImageToRightOriginalNodeMap.find(this_node_index) != mImageToRightOriginalNodeMap.end())
            {
                number_of_right_image_nodes++;
            }
        }

        // Elements on the left hand side (images of right nodes)
        if (number_of_right_image_nodes == 1 || number_of_right_image_nodes == 2)
        {
            mLeftPeriodicBoundaryElementIndices.insert(elem_iter->GetIndex());
        }

        // Elements on the right (images of left nodes)
        if (number_of_left_image_nodes == 1 || number_of_left_image_nodes == 2)
        {
            mRightPeriodicBoundaryElementIndices.insert(elem_iter->GetIndex());
        }
    }
    // Every boundary element on the left should have a corresponding element on
    // the right...
    assert(mLeftPeriodicBoundaryElementIndices.size()==mRightPeriodicBoundaryElementIndices.size());
}


unsigned Cylindrical2dMesh::GetCorrespondingNodeIndex(unsigned nodeIndex)
{
    unsigned corresponding_node_index = UINT_MAX;

    // If nodeIndex is a member of mRightOriginals, then find the corresponding node index in mRightImages
    std::vector<unsigned>::iterator right_orig_iter = std::find(mRightOriginals.begin(), mRightOriginals.end(), nodeIndex);
    if (right_orig_iter != mRightOriginals.end())
    {
        corresponding_node_index = mRightImages[right_orig_iter - mRightOriginals.begin()];
    }
    else
    {
        // If nodeIndex is a member of mRightImages, then find the corresponding node index in mRightOriginals
        std::vector<unsigned>::iterator right_im_iter = std::find(mRightImages.begin(), mRightImages.end(), nodeIndex);
        if (right_im_iter != mRightImages.end())
        {
            corresponding_node_index = mRightOriginals[right_im_iter - mRightImages.begin()];
        }
        else
        {
            // If nodeIndex is a member of mLeftOriginals, then find the corresponding node index in mLeftImages
            std::vector<unsigned>::iterator left_orig_iter = std::find(mLeftOriginals.begin(), mLeftOriginals.end(), nodeIndex);
            if (left_orig_iter != mLeftOriginals.end())
            {
                corresponding_node_index = mLeftImages[left_orig_iter - mLeftOriginals.begin()];
            }
            else
            {
                // If nodeIndex is a member of mLeftImages, then find the corresponding node index in mLeftOriginals
                std::vector<unsigned>::iterator left_im_iter = std::find(mLeftImages.begin(), mLeftImages.end(), nodeIndex);
                if (left_im_iter != mLeftImages.end())
                {
                    corresponding_node_index = mLeftOriginals[left_im_iter - mLeftImages.begin()];
                }
            }
        }
    }

    // We must have found the corresponding node index
    assert(corresponding_node_index != UINT_MAX);
    return corresponding_node_index;
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Cylindrical2dMesh)
