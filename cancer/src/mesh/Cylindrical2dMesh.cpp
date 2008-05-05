#ifndef _CYLINDRICAL2DMESH_CPP_
#define _CYLINDRICAL2DMESH_CPP_

#include "Cylindrical2dMesh.hpp"

/**
 * @param pElement  
 * @param rImageNodes Left or right image nodes
 * @param rOriginalNodes  Left or right original nodes
 * @param nodeIndex 
 */
void Cylindrical2dMesh::ReplaceImageWithRealNodeOnElement(Element<2,2>* pElement, std::vector<unsigned> &rImageNodes, std::vector<unsigned> &rOriginalNodes, unsigned nodeIndex ) 
{ 
    for (unsigned j=0; j<rImageNodes.size(); j++)
    {
        if(nodeIndex==rImageNodes[j])
        {
            pElement->ReplaceNode(mNodes[rImageNodes[j]],mNodes[rOriginalNodes[j]]);
        }
    }
}

/**
 * Constructor
 * 
 * @param width the width of the crypt (circumference) 
 */
Cylindrical2dMesh::Cylindrical2dMesh(double width)
  : ConformingTetrahedralMesh<2, 2>(),
    mWidth(width)    
{
    assert(width > 0.0);
}


Cylindrical2dMesh::Cylindrical2dMesh(double width, std::vector<Node<2> *> nodes)
  : ConformingTetrahedralMesh<2, 2>(),
    mWidth(width)
{
    assert(width > 0.0);
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<2>* temp_node = nodes[index];
        double x = temp_node->rGetLocation()[0];
        x=x; // Fix optimised build
        assert( 0 <= x && x < width);
        mNodes.push_back(temp_node);
    }
    
    NodeMap node_map(nodes.size());
    ReMesh(node_map);
}

/**
 * Calls GetWidthExtremes on the Conforming mesh class
 * to calculate mTop and mBottom for the cylindrical mesh.
 */
void Cylindrical2dMesh::UpdateTopAndBottom()
{
    c_vector<double,2> extremes = GetWidthExtremes(1);
    mBottom = extremes[0];
    mTop = extremes[1];   
}


/**
 * Creates a set of mirror nodes for a cylindrical re-mesh. 
 * All mesh points should be 0<x<mWidth.
 * 
 * Updates mRightImages and mLeftImages
 */
void Cylindrical2dMesh::CreateMirrorNodes()
{
    unsigned num_nodes=GetNumAllNodes();
    double half_way = (mWidth)/2.0;
            
    //TestTopAndBottomRowAlignment();
    mLeftOriginals.clear();
    mLeftImages.clear();
    mRightOriginals.clear();
    mRightImages.clear();
    mLeftPeriodicBoundaryElementIndices.clear();
    mRightPeriodicBoundaryElementIndices.clear();

    for (unsigned i=0; i<num_nodes; i++)
    {
        if (!mNodes[i]->IsDeleted())
        {
            c_vector<double, 2> location = mNodes[i]->rGetLocation();
            unsigned this_node_index = mNodes[i]->GetIndex();
            double this_node_x_location = location[0];
            
            // Check the mesh currently conforms to the dimensions given.        
            assert(0.0<=location[0]);
            assert(location[0]<=mWidth);
            
            // Put the nodes which are to be mirrored in the relevant vectors
            if (this_node_x_location<half_way)
            {
                mLeftOriginals.push_back(this_node_index);
            }
            else
            {
                mRightOriginals.push_back(this_node_index);
            }
        }
    }

    // Go through the left original nodes and create an image node
    // recording its new index.
    for (unsigned i=0; i<mLeftOriginals.size(); i++)
    {
        c_vector<double, 2> location = mNodes[mLeftOriginals[i]]->rGetLocation();
        location[0] = location[0] + mWidth;

        unsigned new_node_index = ConformingTetrahedralMesh<2,2>::AddNode(new Node<2>(0u, location));
        mLeftImages.push_back(new_node_index);
    }
    
    // Go through the right original nodes and create an image node
    // recording its new index.
    for (unsigned i=0; i<mRightOriginals.size(); i++)
    {
        // Create new image nodes
        c_vector<double, 2> location = mNodes[mRightOriginals[i]]->rGetLocation();
        location[0] = location[0] - mWidth;

        unsigned new_node_index = ConformingTetrahedralMesh<2,2>::AddNode(new Node<2>(0u, location));
        mRightImages.push_back(new_node_index);
    }
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
       unsigned new_node_index = ConformingTetrahedralMesh<2,2>::AddNode(new Node<2>(0u, location));
       mTopHaloNodes.push_back(new_node_index);
       
       location[1] = y_bottom_coordinate;
       new_node_index = ConformingTetrahedralMesh<2,2>::AddNode(new Node<2>(0u, location));
       mBottomHaloNodes.push_back(new_node_index);
    }
}

/**
 * Conducts a cylindrical remesh (OVERRIDDEN constructor of main ReMesh function)
 * 
 * Firstly calls CreateMirrorNodes to create mirror image nodes
 * Then calls remesher
 * Maps new node indices
 * calls ReconstructCylindricalMesh to remove surplus nodes to create a fully periodic mesh.
 * 
 * @param &map a reference to a nodemap which should be created with the required number of nodes.
 */
void Cylindrical2dMesh::ReMesh(NodeMap &map)
{
    unsigned old_num_all_nodes = GetNumAllNodes();
	
	map.Resize(old_num_all_nodes);
	map.ResetToIdentity();
    
    // flag the deleted nodes as deleted in the map
	for(unsigned i=0; i<old_num_all_nodes; i++)
	{
		if(mNodes[i]->IsDeleted())
		{
			map.SetDeleted(i);
		}
	}
	
    CreateHaloNodes();
    
    // Create a mirrored load of nodes for the normal remesher to work with.
    CreateMirrorNodes();

    // The mesh now has messed up boundary elements 
    // but this doesn't matter as the ReMesh below
    // doesn't read them in and reconstructs the
    // boundary elements.
    
    // Call the normal re-mesh
    // note that the mesh now has lots of extra nodes which will be deleted, hence the name 'big_map'
    NodeMap big_map(GetNumAllNodes()); 
    ConformingTetrahedralMesh<2,2>::ReMeshWithTriangleLibrary(big_map);
    //ConformingTetrahedralMesh<2,2>::ReMesh(big_map);
    // if the big_map isn't the identity map, the little map ('map') needs to be
    // altered accordingly before being passed to the user. not sure how this all works,
    // so deal with this bridge when we get to it 
    assert(big_map.IsIdentityMap());

    // Re-Index the vectors according to the big nodemap.

    for (unsigned i=0; i<mLeftOriginals.size(); i++)
    {
        mLeftOriginals[i] = big_map.GetNewIndex(mLeftOriginals[i]);
        mLeftImages[i] = big_map.GetNewIndex(mLeftImages[i]);
    }
    
    for (unsigned i=0; i<mRightOriginals.size(); i++)
    {
        mRightOriginals[i] = big_map.GetNewIndex(mRightOriginals[i]);
        mRightImages[i] = big_map.GetNewIndex(mRightImages[i]);
    }
    
    for (unsigned i=0; i<mTopHaloNodes.size(); i++)
    {
        mTopHaloNodes[i] = big_map.GetNewIndex(mTopHaloNodes[i]);
        mBottomHaloNodes[i] = big_map.GetNewIndex(mBottomHaloNodes[i]);
    }
    // This method takes in the double sized mesh, 
    // with its new boundary elements,
    // and removes the relevant nodes, elements and boundary elements
    // to leave a proper periodic mesh.
    
//    TrianglesMeshWriter<2,2> mesh_writer("","halo_problem_mesh");
//    mesh_writer.WriteFilesUsingMesh(*this);
    
    GenerateVectorsOfElementsStraddlingPeriodicBoundaries();
    
    CorrectNonPeriodicMesh();
    
    ReconstructCylindricalMesh();
    
    DeleteHaloNodes();
    
    // Create a random (!) boundary element between two nodes of the first element if it is not deleted.
    // This is a temporary measure to get around reindex crashing when there are no boundary elements ( J. Coopers idea )
    bool boundary_element_made = false;
    unsigned elem_index = 0;
    while (elem_index<GetNumAllElements() && !boundary_element_made)
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
            
        }
        elem_index++;   
    }    

    // now call ReIndex to remove the temporary nodes which are marked as deleted. 
	NodeMap reindex_map(GetNumAllNodes());
    ReIndex(reindex_map);
    assert(!reindex_map.IsIdentityMap());  // maybe don't need this
    
    // go through the reindex map and use it to populate the original NodeMap
    // (the one that is returned to the user)
    for(unsigned i=0; i<map.Size(); i++) // only going up to be size of map, not size of reindex_map
    {
        if(reindex_map.IsDeleted(i))
        {
            // i < num_original_nodes and node is deleted, this should correspond to
            // a node that was labelled as before the remeshing, so should have already
            // been set as deleted in the map above 
            assert(map.IsDeleted(i));
        }
        else
        {
            map.SetNewIndex(i, reindex_map.GetNewIndex(i) );
        }
    }
}


/**
 * Deletes the mirror image nodes, elements and boundary elements
 * created for a cylindrical remesh by cycling through the 
 * elements and changing elements with
 * partly real and partly imaginary elements to be real with
 * periodic real nodes instead of mirror image nodes.
 * 
 */
void Cylindrical2dMesh::ReconstructCylindricalMesh()
{
//    TrianglesMeshWriter<2,2> mesh_writer("","debug_periodic_mesh");
//    mesh_writer.WriteFilesUsingMesh(*this);

    // Figure out which elements have real nodes and image nodes in them
    // and replace image nodes with corresponding real ones.
    for (unsigned elem_index = 0; elem_index<GetNumAllElements(); elem_index++)
    {
        Element<2,2>* p_element = GetElement(elem_index);
        if (!p_element->IsDeleted())
        {
            // left images are on the right of the mesh
            unsigned number_of_left_image_nodes = 0u;
            unsigned number_of_right_image_nodes = 0u;
            for (unsigned i=0; i<3; i++)
            {
                unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                //std::cout << "Node " << this_node_index << "\t";
                bool this_node_a_left_image = false;
                bool this_node_a_right_image = false;
                if(IsThisIndexInList(this_node_index,mLeftImages))
                {
                    this_node_a_left_image = true;
//                    std::cout << "Node " << this_node_index << " is a left image\n";
                }
                if(IsThisIndexInList(this_node_index,mRightImages))
                {
                    this_node_a_right_image = true;
//                    std::cout << "Node " << this_node_index << " is a right image\n";
                }
                if(this_node_a_left_image)
                {
                    number_of_left_image_nodes++;
                }
                if(this_node_a_right_image)
                {
                    number_of_right_image_nodes++;
                }
            }
            //std::cout << "Element found to be " << p_element->GetIndex() << " num images of left " << number_of_left_image_nodes <<" num images of right " << number_of_right_image_nodes << "\n" << std::flush;
            // Delete all the elements on the left hand side (images of right)...
            if (number_of_right_image_nodes>=1u)
            {
                //std::cout << "\tright image element (on left) ... deleting\n" << std::flush;
                p_element->MarkAsDeleted();
                mDeletedElementIndices.push_back(p_element->GetIndex());
            }
            
            // Delete only purely imaginary elements on the right (images of left nodes)
            if (number_of_left_image_nodes==3u)
            {
                //std::cout << "\tpurely image element (on right) ... deleting\n" << std::flush;
                p_element->MarkAsDeleted();
                mDeletedElementIndices.push_back(p_element->GetIndex());
            }
            
            /* 
             * If some are images then replace them with the real nodes.
             * 
             * There can be elements with either two image nodes on the right (and one real)
             * or one image node on the right (and two real).
             */
            if (number_of_left_image_nodes==1u || number_of_left_image_nodes==2u )
            {   //
                //std::cout << "\tPeriodic element with nodes\n" << std::flush;   
                for (unsigned i=0; i<3; i++)
                {
                    unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                    //std::cout << "this node index " << this_node_index << " replaced with ";
                    ReplaceImageWithRealNodeOnElement(p_element,mLeftImages,mLeftOriginals,this_node_index);
                    //std::cout << p_element->GetNodeGlobalIndex(i) << "\n" << std::flush;
                }
            }
        }
    }// end of loop over elements

    // Figure out which boundary elements have real nodes and image nodes in them
    // and replace image nodes with corresponding real ones.
    for (unsigned elem_index = 0; elem_index<GetNumAllBoundaryElements(); elem_index++)
    {
        BoundaryElement<1,2>* p_boundary_element = GetBoundaryElement(elem_index);
        if (!p_boundary_element->IsDeleted())
        {
            unsigned number_of_image_nodes = 0;
            for (unsigned i=0; i<2; i++)
            {
                unsigned this_node_index = p_boundary_element->GetNodeGlobalIndex(i);
                //std::cout << this_node_index << "\t";
                bool this_node_an_image = false;
                
                if(IsThisIndexInList(this_node_index,mLeftImages))
                {
                    this_node_an_image = true;
                }
                if(IsThisIndexInList(this_node_index,mRightImages))
                {
                    this_node_an_image = true;
                }
                if(this_node_an_image)
                {
                    number_of_image_nodes++;
                }
            }
                        
            if (number_of_image_nodes==2 )
            {
                //std::cout << "IMAGE\n" << std::flush;
                p_boundary_element->MarkAsDeleted();
                mDeletedBoundaryElementIndices.push_back(p_boundary_element->GetIndex());
            }
            /*
             * To avoid having two copies of the boundary elements on the periodic
             * boundaries we only deal with the elements on the left image and 
             * delete the ones on the right image.
             */
            if (number_of_image_nodes==1 )
            {
                 
                for (unsigned i=0; i<2; i++)
                {
                    unsigned this_node_index = p_boundary_element->GetNodeGlobalIndex(i);
                    for (unsigned j=0; j<mLeftImages.size(); j++)
                    {
                        if(this_node_index==mLeftImages[j])
                        {
                            //std::cout << "PERIODIC \n" << std::flush;  
                            p_boundary_element->ReplaceNode(mNodes[mLeftImages[j]],mNodes[mLeftOriginals[j]]);
                            //std::cout << "Node " << mLeftImages[j] << " swapped for node " << mLeftOriginals[j] << "\n" << std::flush;
                        }
                    }
                    for (unsigned j=0; j<mRightImages.size(); j++)
                    {
                        if(this_node_index==mRightImages[j])
                        {
                            //std::cout << "IMAGE\n" << std::flush;
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


/**
 * This OVERRIDDEN method evaluates the (surface) distance between two points in a 2D Cylindrical geometry.
 * 
 * @param rLocation1 the x and y co-ordinates of point 1
 * @param rLocation2 the x and y co-ordinates of point 2
 * 
 * @return the vector from location1 to location2
 */
c_vector<double, 2> Cylindrical2dMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth>0.0);
    
    c_vector<double, 2> location1 = rLocation1;
    c_vector<double, 2> location2 = rLocation2;
    
    location1[0] = fmod(location1[0], mWidth);
    location2[0] = fmod(location2[0], mWidth);
    
    c_vector<double, 2> vector = location2 - location1;
            
    // handle the cylindrical condition here
    // if the points are more than halfway around the cylinder apart
    // measure the other way.
    if ( vector[0] > (mWidth / 2.0) )
    {
        vector[0] -= mWidth;
    }
    if ( vector[0] < -(mWidth / 2.0))
    {
        vector[0] += mWidth;  
    }
    return vector;
}

/**
 * OVERRIDDEN function to set the location of a node.
 * 
 * If the location should be set outside a cylindrical boundary
 * move it back onto the cylinder.
 * 
 * SetNode moves the node with a particular index to a new point in space and
 * verifies that the signed areas of the supporting Elements are positive
 * @param index is the index of the node to be moved
 * @param point is the new target location of the node
 * @param concreteMove is set to false if we want to skip the signed area tests
 *
 */
void Cylindrical2dMesh::SetNode(unsigned index, ChastePoint<2> point, bool concreteMove)
{
    // Perform a periodic movement if necessary
    if (point.rGetLocation()[0] >= mWidth)
    {   // move point to the left
        point.SetCoordinate(0u, point.rGetLocation()[0]-mWidth);
        //std::cout << "Moving point to the left\n" << std::flush;
    }
    if (point.rGetLocation()[0] < 0.0)
    {   // move point to the right
        point.SetCoordinate(0u, point.rGetLocation()[0]+mWidth);
        //std::cout << "Moving point to the right\n" << std::flush;
    }
    
    // Update the node's location
    ConformingTetrahedralMesh<2,2>::SetNode(index, point, concreteMove); 
}

/**
 * Returns true if an unsigned is contained in a vector of unsigneds
 * 
 * @param rNodeIndex an unsigned value
 * @param rListOfNodes a list of unsigned values
 * 
 * @return whether the unsigned is in this std::vector
 */
bool Cylindrical2dMesh::IsThisIndexInList(const unsigned& rNodeIndex, const std::vector<unsigned>& rListOfNodes)
{
    for (unsigned i=0; i<rListOfNodes.size(); i++)
    {
        if(rNodeIndex==rListOfNodes[i])
        {
            return true;
        }
    }
    return false;
}

/**
 * OVERRIDDEN FUNCTION
 * @param rDimension must be 0 (x) or 1 (y)
 * @return width the CryptWidth or current height 
 */
double Cylindrical2dMesh::GetWidth(const unsigned& rDimension) const
{
    double width=0.0;
    assert(rDimension==0 || rDimension==1);
    if (rDimension==0)
    {
        width = mWidth;   
    }
    else
    {
        width = ConformingTetrahedralMesh<2,2>::GetWidth(rDimension);
    }
    return width;   
}


/** 
 * Add a node to the mesh.
 * 
 * NB. After calling this one or more times, you must then call ReMesh
 *
 */
unsigned Cylindrical2dMesh::AddNode(Node<2> *pNewNode)
{
    unsigned node_index = ConformingTetrahedralMesh<2,2>::AddNode(pNewNode);
    
    // If necessary move it to be back on the cylinder
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point, false); 
    
    return node_index;
}

void Cylindrical2dMesh::CorrectNonPeriodicMesh()
{
//    TrianglesMeshWriter<2,2> mesh_writer("","debug_periodic_mesh");
//    mesh_writer.WriteFilesUsingMesh(*this);
    
    /* 
     * Copy the member variables into new vectors which we modify by knocking out 
     * elements which pair up on each side
     */
    std::set<unsigned> temp_left_hand_side_elements = mLeftPeriodicBoundaryElementIndices;
    std::set<unsigned> temp_right_hand_side_elements = mRightPeriodicBoundaryElementIndices;
    
    for (std::set<unsigned>::iterator left_iter = mLeftPeriodicBoundaryElementIndices.begin(); 
         left_iter != mLeftPeriodicBoundaryElementIndices.end(); 
         ++left_iter)
    {
        unsigned elem_index = *left_iter;
        Element<2,2>* p_element = GetElement(elem_index);
        
        c_vector<unsigned,3> original_element_node_indices;
        c_vector<unsigned,3> corresponding_element_node_indices;
        for (unsigned i=0 ; i<3 ; i++)
        {
            original_element_node_indices[i] = p_element->GetNodeGlobalIndex(i);
            corresponding_element_node_indices[i] = GetCorrespondingNodeIndex(original_element_node_indices[i]);
        }
    
        // search the right hand sides for the coresponding element 
        for (std::set<unsigned>::iterator right_iter = mRightPeriodicBoundaryElementIndices.begin(); 
             right_iter != mRightPeriodicBoundaryElementIndices.end(); 
             ++right_iter)
        {
            unsigned corresponding_elem_index = *right_iter;
            Element<2,2>* p_corresponding_element = GetElement(corresponding_elem_index);
            
            bool is_coresponding_node = true;
                 
            for( unsigned i=0; i<3; i++)
            {
                if( !(corresponding_element_node_indices[i] == p_corresponding_element->GetNodeGlobalIndex(0)) &&
                    !(corresponding_element_node_indices[i] == p_corresponding_element->GetNodeGlobalIndex(1)) &&
                    !(corresponding_element_node_indices[i] == p_corresponding_element->GetNodeGlobalIndex(2)) )
                {
                    is_coresponding_node=false;
                }
            }
            
            if (is_coresponding_node)
            {
                // remove original and coresponding element from sets
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
    assert(temp_left_hand_side_elements.size()<=2u);   
    assert(temp_right_hand_side_elements.size()<=2u);  

    /*
     * Now we just have to use the first pair of elements and copy their info over to the other side.
     * 
     * First we need to get hold of both elements on either side.
     */
    if (temp_left_hand_side_elements.size()==0u || temp_right_hand_side_elements.size()==0u)
    {   
        assert(temp_right_hand_side_elements.size()==0u);
        assert(temp_left_hand_side_elements.size()==0u);
        //std::cout << "No problem elements\n" << std::flush;
    }
    else
    {
        //std::cout << "Problem elements\n" << std::flush;
        if (temp_right_hand_side_elements.size()==2u)
        {   // Use the right hand side meshing and map to left
            #define COVERAGE_IGNORE
            UseTheseElementsToDecideMeshing(temp_right_hand_side_elements);
            #undef COVERAGE_IGNORE
            
        }
        else if (temp_left_hand_side_elements.size()==2u)
        {   // Use the left hand side meshing and map to right
            UseTheseElementsToDecideMeshing(temp_left_hand_side_elements);
        }
        else
        {   // If you get here there are more than two mixed up elements on the periodic edge.
            NEVER_REACHED;
        }
    }
}

void Cylindrical2dMesh::UseTheseElementsToDecideMeshing(std::set<unsigned> mainSideElements)
{
    assert(mainSideElements.size()==2u);
    
//    TrianglesMeshWriter<2,2> mesh_writer("","debug_periodic_mesh");
//    mesh_writer.WriteFilesUsingMesh(*this);
    
    // We find the four nodes surrounding the dodgy meshing, on each side.
    std::set<unsigned> main_four_nodes;
    for (std::set<unsigned>::iterator left_iter = mainSideElements.begin(); 
         left_iter != mainSideElements.end(); 
         ++left_iter)
    {
        unsigned elem_index = *left_iter;
        Element<2,2>* p_element = GetElement(elem_index);
        for (unsigned i=0; i<3 ; i++)
        {
            unsigned index = p_element->GetNodeGlobalIndex(i);
            main_four_nodes.insert(index);
        }
    }
    assert(main_four_nodes.size()==4u);
    
    std::set<unsigned> other_four_nodes;
    for (std::set<unsigned>::iterator iter = main_four_nodes.begin(); 
         iter != main_four_nodes.end(); 
         ++iter)
    {
        other_four_nodes.insert(GetCorrespondingNodeIndex(*iter));
    }
    assert(other_four_nodes.size()==4u);
    
    
    // Find the elements surrounded by the nodes on the right 
    // and change them to match the elements on the left
    std::vector<unsigned> corresponding_elements;
    // loop over all elements
    for (unsigned elem_index = 0; elem_index<GetNumAllElements(); elem_index++)
    {
        Element<2,2>* p_element = GetElement(elem_index);
        if (!p_element->IsDeleted())
        {
        // loop over the nodes of the element
            if (!(other_four_nodes.find(p_element->GetNodeGlobalIndex(0))==other_four_nodes.end()) &&
                !(other_four_nodes.find(p_element->GetNodeGlobalIndex(1))==other_four_nodes.end()) &&
                !(other_four_nodes.find(p_element->GetNodeGlobalIndex(2))==other_four_nodes.end()) )
            {
                corresponding_elements.push_back(elem_index);
                p_element->MarkAsDeleted();
                mDeletedElementIndices.push_back(p_element->GetIndex());
            }
        }
    }
    assert(corresponding_elements.size()==2u);
    
    // Now corresponding_elements contains the two elements which are going to be replaced by mainSideElements
    for (std::set<unsigned>::iterator iter = mainSideElements.begin(); 
         iter != mainSideElements.end(); 
         ++iter)
    {
        Element<2,2>* p_main_element = GetElement(*iter);
        std::vector<Node<2>*> nodes;
        // Put corresponding nodes into a std::vector.
        for (unsigned i=0; i<3 ; i++)
        {
            unsigned main_node = p_main_element->GetNodeGlobalIndex(i);
            nodes.push_back(this->GetNode(GetCorrespondingNodeIndex(main_node)));
        }
        
        // Make a new element.                
        Element<2,2>* p_new_element = new Element<2,2>(GetNumAllElements(), nodes);
        this->mElements.push_back(p_new_element);
    }
    
    // Reindex to get rid of extra elements indices...
    NodeMap map(GetNumAllNodes());
    this->ReIndex(map);    
    
//    TrianglesMeshWriter<2,2> mesh_writer2("","debug_periodic_mesh.1");
//    mesh_writer2.WriteFilesUsingMesh(*this);
}

void Cylindrical2dMesh::GenerateVectorsOfElementsStraddlingPeriodicBoundaries()
{
    mLeftPeriodicBoundaryElementIndices.clear();
    mRightPeriodicBoundaryElementIndices.clear();
    
    for (unsigned elem_index = 0; elem_index<GetNumAllElements(); elem_index++)
    {
        Element<2,2>* p_element = GetElement(elem_index);
        if (!p_element->IsDeleted())
        {
            // left images are on the right of the mesh
            unsigned number_of_left_image_nodes = 0u;
            unsigned number_of_right_image_nodes = 0u;
            for (unsigned i=0; i<3; i++)
            {
                unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                //std::cout << "Node " << this_node_index << "\t";
                bool this_node_a_left_image = false;
                bool this_node_a_right_image = false;
                if(IsThisIndexInList(this_node_index,mLeftImages))
                {
                    this_node_a_left_image = true;
                }
                if(IsThisIndexInList(this_node_index,mRightImages))
                {
                    this_node_a_right_image = true;
                }
                if(this_node_a_left_image)
                {
                    number_of_left_image_nodes++;
                }
                if(this_node_a_right_image)
                {
                    number_of_right_image_nodes++;
                }
            }
            // elements on the left hand side (images of right)...
            if (number_of_right_image_nodes==1u || number_of_right_image_nodes==2u)
            {
                mLeftPeriodicBoundaryElementIndices.insert(elem_index);
            }
            
            // elements on the right (images of left nodes)
            if (number_of_left_image_nodes==1u|| number_of_left_image_nodes==2u)
            {
                mRightPeriodicBoundaryElementIndices.insert(elem_index);
            }
        }
    }
}

unsigned Cylindrical2dMesh::GetCorrespondingNodeIndex(unsigned nodeIndex)
{
    unsigned corresponding_node_index = UINT_MAX;
    bool found = false;
    
    if (IsThisIndexInList(nodeIndex, mRightOriginals))
    {
        for (unsigned i=0 ; i<mRightOriginals.size() ; i++)
        {
            if (mRightOriginals[i]==nodeIndex)
            {
                corresponding_node_index = mRightImages[i];
                found = true;
            }   
        }   
    }
    if (IsThisIndexInList(nodeIndex, mRightImages))
    {
        for (unsigned i=0 ; i<mRightImages.size() ; i++)
        {
            if (mRightImages[i]==nodeIndex)
            {
                corresponding_node_index = mRightOriginals[i];
                found = true;
            }   
        }   
    }
    
    if (IsThisIndexInList(nodeIndex, mLeftOriginals))
    {
        for (unsigned i=0 ; i<mLeftOriginals.size() ; i++)
        {
            if (mLeftOriginals[i]==nodeIndex)
            {
                corresponding_node_index = mLeftImages[i];
                found = true;
            }   
        }   
    }    
    if (IsThisIndexInList(nodeIndex, mLeftImages))
    {
        for (unsigned i=0 ; i<mLeftImages.size() ; i++)
        {
            if (mLeftImages[i]==nodeIndex)
            {
                corresponding_node_index = mLeftOriginals[i];
                found = true;
            }   
        }   
    }
    
    assert(found);
    return corresponding_node_index;
}

#endif //_CYLINDRICAL2DMESH_CPP_
