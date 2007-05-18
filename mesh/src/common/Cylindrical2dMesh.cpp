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
    for (unsigned j=0 ; j<rImageNodes.size() ; j++)
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
 * @param xTop the y-co-ord of the top row of ghost nodes (important to keep them together)
 * @param xBottom the y-coord of the bottom row of ghost nodes
 * @param topBoundary the inidces of the top ghost nodes
 * @param bottomBoundary the inidces of the bottom ghost nodes
 */
Cylindrical2dMesh::Cylindrical2dMesh(double width, double xTop, double xBottom, std::vector<unsigned > topBoundary, std::vector<unsigned > bottomBoundary)
  : ConformingTetrahedralMesh<2, 2>(),
    mWidth(width), 
    mTop(xTop),
    mBottom(xBottom),
    mTopBoundary(topBoundary),
    mBottomBoundary(bottomBoundary)
{
    
    assert(width > 0.0);
    
}
    


/**
 * Creates a set of mirror nodes for a cylindrical re-mesh. 
 * All mesh points should be 0<x<mWidth.
 * 
 * Updates mRightImages and mLeftImages
 */
void Cylindrical2dMesh::CreateMirrorNodes()
{
    unsigned num_nodes=GetNumNodes();
    double half_way = (mWidth)/2.0;
            
    TestTopAndBottomRowAlignment();
    
    mLeftOriginals.clear();
    mLeftImages.clear();
    mRightOriginals.clear();
    mRightImages.clear();
    
    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, 2> location = mNodes[i]->rGetLocation();
        unsigned this_node_index = mNodes[i]->GetIndex();
        double this_node_x_location = location[0];
        
        // Check the mesh currently conforms to the dimensions given.        
        assert(0.0<=location[0] && location[0]<=mWidth);
        
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
    
    // Go through the left original nodes and create an image node
    // recording its new index.
    for (unsigned i=0 ; i<mLeftOriginals.size() ; i++)
    {
        c_vector<double, 2> location = mNodes[mLeftOriginals[i]]->rGetLocation();
        location[0] = location[0] + mWidth;

        unsigned new_node_index = ConformingTetrahedralMesh<2,2>::AddNode(new Node<2>(0u, location));
        mLeftImages.push_back(new_node_index);
    }
    
    // Go through the right original nodes and create an image node
    // recording its new index.
    for (unsigned i=0 ; i<mRightOriginals.size() ; i++)
    {
        // Create new image nodes
        c_vector<double, 2> location = mNodes[mRightOriginals[i]]->rGetLocation();
        location[0] = location[0] - mWidth;

        unsigned new_node_index = ConformingTetrahedralMesh<2,2>::AddNode(new Node<2>(0u, location));
        mRightImages.push_back(new_node_index);
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
    // Create a mirrored load of nodes for the normal remesher to work with.
    CreateMirrorNodes();

    // The mesh now has messed up boundary elements 
    // but this doesn't matter as the ReMesh below
    // doesn't read them in and reconstructs the
    // boundary elements.
    
    // Call the normal re-mesh
    ConformingTetrahedralMesh<2,2>::ReMesh(map);
    
    //
    // Re-Index the vectors according to the node_map.
    //
    for (unsigned i = 0 ; i<mLeftOriginals.size() ; i++)
    {
            mLeftOriginals[i]=map.GetNewIndex(mLeftOriginals[i]);
            mLeftImages[i]=map.GetNewIndex(mLeftImages[i]);
    }
    for (unsigned i = 0 ; i<mRightOriginals.size() ; i++)
    {
            mRightOriginals[i]=map.GetNewIndex(mRightOriginals[i]);
            mRightImages[i]=map.GetNewIndex(mRightImages[i]);
    }
    for (unsigned i = 0 ; i<mTopBoundary.size() ; i++)
    {
            mTopBoundary[i]=map.GetNewIndex(mTopBoundary[i]);
    }
    for (unsigned i = 0 ; i<mBottomBoundary.size() ; i++)
    {
            mBottomBoundary[i]=map.GetNewIndex(mBottomBoundary[i]);
    }
    
    // This method takes in the double sized mesh, 
    // with its new boundary elements,
    // and removes the relevant nodes, elements and boundary elements
    // to leave a proper periodic mesh.
    ReconstructCylindricalMesh();
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
    
    // Figure out which elements have real nodes and image nodes in them
    // and replace image nodes with corresponding real ones.
    for (unsigned elem_index = 0; elem_index<GetNumAllElements(); elem_index++)
    {
        Element<2,2>* p_element = GetElement(elem_index);
        if (!p_element->IsDeleted())
        {
            unsigned number_of_image_nodes = 0;
            for (unsigned i=0 ; i<3 ; i++)
            {
                unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                //std::cout << "Node " << this_node_index << "\t";
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
            
            //std::cout << "\nNumber of image nodes = " << number_of_image_nodes << "\n" << std::flush;
            if (number_of_image_nodes==3 || number_of_image_nodes==2)
            {
                //std::cout << "purely image element\n" << std::flush;
                p_element->MarkAsDeleted();
                mDeletedElementIndices.push_back(p_element->GetIndex());
            }
            /* 
             * If some are images then replace them with the real nodes.
             * 
             * There would be two copies of each periodic element
             * one with one image and two real (on one side),
             * another with two images and one real (on the other side).
             * Because of this we can just take one case (one image node)
             * and delete the other element
             */
            if (number_of_image_nodes==1)
            {
                
                //std::cout << "Periodic element found \n" << std::flush;   
                for (unsigned i=0 ; i<3 ; i++)
                {
                    unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                    ReplaceImageWithRealNodeOnElement(p_element,mLeftImages,mLeftOriginals,this_node_index);
                    ReplaceImageWithRealNodeOnElement(p_element,mRightImages,mRightOriginals,this_node_index);
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
            //std::cout << "Boundary Element " << elem_index << " connects nodes : ";
            unsigned number_of_image_nodes = 0;
            for (unsigned i=0 ; i<2 ; i++)
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
                        
            if (number_of_image_nodes==2)
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
            if (number_of_image_nodes==1)
            {
                 
                for (unsigned i=0 ; i<2 ; i++)
                {
                    unsigned this_node_index = p_boundary_element->GetNodeGlobalIndex(i);
                    for (unsigned j=0 ; j<mLeftImages.size() ; j++)
                    {
                        if(this_node_index==mLeftImages[j])
                        {
                            //std::cout << "PERIODIC \n" << std::flush;  
                            p_boundary_element->ReplaceNode(mNodes[mLeftImages[j]],mNodes[mLeftOriginals[j]]);
                            //std::cout << "Node " << mLeftImages[j] << " swapped for node " << mLeftOriginals[j] << "\n" << std::flush;
                        }
                    }
                    for (unsigned j=0 ; j<mRightImages.size() ; j++)
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
    
    // Delete all image nodes
    for (unsigned i=0 ; i<mLeftImages.size() ; i++)
    {
        mNodes[mLeftImages[i]]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(mLeftImages[i]);
    }
    for (unsigned i=0 ; i<mRightImages.size() ; i++)
    {
        mNodes[mRightImages[i]]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(mRightImages[i]);
    }
            
    // ReIndex the mesh
    ReIndex();
}

/**
 * This OVERRIDDEN method evaluates the (surface) distance between two points in a 2D Cylindrical geometry.
 * 
 * locations should lie between [0, mWidth) 
 * 
 * @param rLocation1 the x and y co-ordinates of point 1
 * @param rLocation2 the x and y co-ordinates of point 2
 * 
 * @return the vector from location1 to location2
 */
c_vector<double, 2> Cylindrical2dMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth>0.0);
    assert(0.0<=rLocation1[0]);  // 1st point is not in cylinder
    assert(0.0<=rLocation2[0]);  // 2nd point is not in cylinder
    assert(mWidth>=rLocation1[0]);  // 1st point is not in cylinder
    assert(mWidth>=rLocation2[0]);  // 2nd point is not in cylinder
    
    c_vector<double, 2> vector = rLocation2 - rLocation1;
            
    // handle the cylindrical condition here
    // if the points are more than halfway around the cylinder apart
    // measure the other way.
    if ( vector(0) > (mWidth / 2.0) )
    {
        vector(0) -= mWidth;
    }
    if ( vector(0) < -(mWidth / 2.0))
    {
        vector(0) += mWidth;  
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
void Cylindrical2dMesh::SetNode(unsigned index, Point<2> point, bool concreteMove)
{
    // We need to move all of the nodes in the top and bottom boundaries together.
    bool on_the_top_of_cylinder = IsThisIndexInList(index,mTopBoundary);
    bool on_the_bottom_of_cylinder = IsThisIndexInList(index,mBottomBoundary);
    
    if (on_the_top_of_cylinder || on_the_bottom_of_cylinder)
    {
        double y_co_ord = point.rGetLocation()[1];
        if(on_the_top_of_cylinder)
        {
            for (unsigned i=0 ; i<mTopBoundary.size() ; i++)
            {
                // Get each node's x position and update so that y position matches on each
                Point<2> boundary_point = ConformingTetrahedralMesh<2,2>::mNodes[mTopBoundary[i]]->GetPoint();
                boundary_point.SetCoordinate(1u, y_co_ord);
                ConformingTetrahedralMesh<2,2>::SetNode(mTopBoundary[i], boundary_point, concreteMove); 
            }
            mTop = y_co_ord; // update the boundary definition.
        }
        if(on_the_bottom_of_cylinder)
        {
            for (unsigned i=0 ; i<mBottomBoundary.size() ; i++)
            {
                // Get each node's x position and update so that y position matches on each
                Point<2> boundary_point = ConformingTetrahedralMesh<2,2>::mNodes[mBottomBoundary[i]]->GetPoint();
                boundary_point.SetCoordinate(1u, y_co_ord);
                ConformingTetrahedralMesh<2,2>::SetNode(mBottomBoundary[i], boundary_point, concreteMove); 
            }
            mBottom = y_co_ord; // update the boundary definition.
        }
        
    }
    
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
    bool is_in_vector = false;
    for (unsigned i=0 ; i<rListOfNodes.size() ; i++)
    {
        if(rNodeIndex==rListOfNodes[i])
        {
            is_in_vector = true;
        }
    }
    return is_in_vector;
}

void Cylindrical2dMesh::TestTopAndBottomRowAlignment()
{
    // Check that the top and bottom rows have the same y-co-ordinate 
    // or things will start to go wrong with boundary elements.
    double y_location = 0.0;
    for (unsigned i=0 ; i<mTopBoundary.size() ; i++)
    {
        y_location = mNodes[mTopBoundary[i]]->rGetLocation()[1];
        if (fabs(y_location - mTop)>1e-3)
        {
            std::cout << "y = " << y_location << ", mTop = " << mTop << "\n" << std::flush;
            EXCEPTION("The top row of ghost nodes is not aligned.");   
        }
    }
    for (unsigned i=0 ; i<mBottomBoundary.size() ; i++)
    {
        y_location = mNodes[mBottomBoundary[i]]->rGetLocation()[1];
        if (fabs(y_location - mBottom)>1e-3)
        {
            EXCEPTION("The bottom row of ghost nodes is not aligned.");   
        }
    }
}

/**
 * OVERRIDDEN FUNCTION
 * @param rDimension must be 0 (x) or 1 (y)
 * @return width the CryptWidth or current height 
 */
double Cylindrical2dMesh::GetWidth(const unsigned& rDimension)
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
    Point<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point, false); 
    
    return node_index;
}
    

#endif //_CYLINDRICAL2DMESH_CPP_
