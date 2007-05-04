#ifndef _CYLINDRICAL2DMESH_CPP_
#define _CYLINDRICAL2DMESH_CPP_

#include "ConformingTetrahedralMesh.cpp"
#include "NodeMap.hpp"

class Cylindrical2dMesh : public ConformingTetrahedralMesh<2, 2>
{
private:
    /** The left x-coord of cylinder boundary*/
    double mXLeft;
    /** The right x-coord of the same cylinder boundary*/
    double mXRight;
    /** The circumference of the cylinder */
    double mWidth;
    /** The top of the cylinder (y-coord) */
    double mTop;
    /** The bottom of the cylinder (y-coord) */
    double mBottom;
    
    /** The indices of nodes on the top boundary */
    std::vector<unsigned > mTopBoundary;
    /** The indices of nodes on the bottom boundary */
    std::vector<unsigned > mBottomBoundary;
    
public:
    
    /**
     * Constructor
     * 
     * @param x0 the left hand cylindrical boundary (usually 0.0)
     * @param x1 the right hand cylindrical boundary (usually CryptWidth)
     * @param xTop the y-co-ord of the top row of ghost nodes (important to keep them together)
     * @param xBottom the y-coord of the bottom row of ghost nodes
     * @param topBoundary the inidces of the top ghost nodes
     * @param bottomBoundary the inidces of the bottom ghost nodes
     */
    Cylindrical2dMesh(double x0, double x1, double xTop, double xBottom, std::vector<unsigned > topBoundary, std::vector<unsigned > bottomBoundary): ConformingTetrahedralMesh<2, 2>()
    {
        assert(x1>x0);
        mXLeft = x0;
        mXRight = x1;
        mWidth = mXRight - mXLeft;
        mTop = xTop;
        mBottom = xBottom;
        mTopBoundary = topBoundary;
        mBottomBoundary = bottomBoundary;
    }
    
    Cylindrical2dMesh(): ConformingTetrahedralMesh<2, 2>()
    {
        EXCEPTION("Please specify the boundaries of the cylinder using the other constructor");
    }
    
    ~Cylindrical2dMesh()
    {
        
    }
    
    /**
     * Creates a set of mirror nodes for a cylindrical re-mesh. 
     * Where x1 - x0 = 2*PI in a cylinder. 
     * All mesh points should be x0<x<x1.
     * 
     * @return a map of four standard vectors of indices
     * 0. The left nodes which have been mirrored
     * 1. The image nodes relating to these left nodes (on right of mesh)
     * 2. The right nodes which have been mirrored
     * 3. The image nodes relating to these right nodes (on left of mesh)
     */
    std::vector<std::vector <unsigned> > CreateMirrorNodes()
    {
        unsigned num_nodes=GetNumNodes();
        double half_way = (mWidth)/2.0;
                
        TestTopAndBottomRowAlignment();
        
        // Label the nodes which will have to be mirrored
        std::vector<unsigned> left_original_node_indices;
        std::vector<unsigned> lefts_image_node_indices;
        std::vector<unsigned> right_original_node_indices;
        std::vector<unsigned> rights_image_node_indices;
     
        for (unsigned i=0; i<num_nodes; i++)
        {
            c_vector<double, 2> location = mNodes[i]->rGetLocation();
            unsigned this_node_index = mNodes[i]->GetIndex();
            double this_node_x_location = location[0];
            
            // Check the mesh currently conforms to the dimensions given.
            if (!(mXLeft<=location[0] && location[0]<=mXRight))
            {
#define COVERAGE_IGNORE
                std::cout << "node(" << i << "), x = " << location[0] << "\n" << std::flush;
                std::cout << "left boundary = " << mXLeft << ", right boundary = " << mXRight << std::flush;
                EXCEPTION("A node lies outside the cylindrical region");    
#undef COVERAGE_IGNORE
            }
            
            // Put the nodes which are to be mirrored in the relevant vectors
            if (this_node_x_location<half_way)
            {
                left_original_node_indices.push_back(this_node_index);
            }
            if (this_node_x_location>=half_way)
            {
                right_original_node_indices.push_back(this_node_index);
            }
        }
        
        // Go through the left original nodes and create an image node
        // recording its new index.
        for (unsigned i=0 ; i<left_original_node_indices.size() ; i++)
        {
            c_vector<double, 2> location = mNodes[left_original_node_indices[i]]->rGetLocation();
            location[0] = location[0] + mWidth;
    
            unsigned new_node_index = AddNode(new Node<2>(0u, location));
            lefts_image_node_indices.push_back(new_node_index);
        }
        
        // Go through the right original nodes and create an image node
        // recording its new index.
        for (unsigned i=0 ; i<right_original_node_indices.size() ; i++)
        {
            // Create new image nodes
            c_vector<double, 2> location = mNodes[right_original_node_indices[i]]->rGetLocation();
            location[0] = location[0] - mWidth;
    
            unsigned new_node_index = AddNode(new Node<2>(0u, location));
            rights_image_node_indices.push_back(new_node_index);
        }
        
        // Return all the information...
        std::vector<std::vector <unsigned> > image_map;
        image_map.push_back(left_original_node_indices);
        image_map.push_back(lefts_image_node_indices);
        image_map.push_back(right_original_node_indices);
        image_map.push_back(rights_image_node_indices);
        return image_map;
    }

    /**
     * Conducts a cylindrical remesh (overwritten constructor of main ReMesh function)
     * 
     * Firstly calls CreateMirrorNodes to create mirror image nodes
     * Then calls remesher
     * Maps new node indices
     * calls ReconstructCylindricalMesh to remove surplus nodes to create a fully periodic mesh.
     * 
     * @param &map a reference to a nodemap which should be created with the required number of nodes.
     */
    void ReMesh(NodeMap &map)
    {
        // std::cout << "Conducting a cylindrical re-mesh\n" << std::flush;
        // Create a mirrored load of nodes for the normal remesher to work with.
        std::vector<std::vector<unsigned> > image_map = CreateMirrorNodes();
    
        // The mesh now has messed up boundary elements 
        // but this doesn't matter as the ReMesh below
        // doesn't read them in and reconstructs the
        // boundary elements.
    
        std::vector<unsigned> left_original = image_map[0];
        std::vector<unsigned> left_images = image_map[1];
        std::vector<unsigned> right_original = image_map[2];
        std::vector<unsigned> right_images = image_map[3];
        
        // Call the normal re-mesh
        ConformingTetrahedralMesh<2,2>::ReMesh(map);
        
        //
        // Re-Index the vectors according to the node_map.
        //
        for (unsigned i = 0 ; i<left_original.size() ; i++)
        {
                left_original[i]=map.GetNewIndex(left_original[i]);
                left_images[i]=map.GetNewIndex(left_images[i]);
        }
        for (unsigned i = 0 ; i<right_original.size() ; i++)
        {
                right_original[i]=map.GetNewIndex(right_original[i]);
                right_images[i]=map.GetNewIndex(right_images[i]);
        }
        for (unsigned i = 0 ; i<mTopBoundary.size() ; i++)
        {
                mTopBoundary[i]=map.GetNewIndex(mTopBoundary[i]);
        }
        for (unsigned i = 0 ; i<mBottomBoundary.size() ; i++)
        {
                mBottomBoundary[i]=map.GetNewIndex(mBottomBoundary[i]);
        }
        
        image_map[0] = left_original;
        image_map[1] = left_images;
        image_map[2] = right_original;
        image_map[3] = right_images;
        
        // This method takes in the double sized mesh, 
        // with its new boundary elements,
        // and removes the relevant nodes, elements and boundary elements
        // to leave a proper periodic mesh.
        ReconstructCylindricalMesh(image_map);
    }

    /**
     * Deletes the mirror image nodes, elements and boundary elements
     * created for a cylindrical remesh by cycling through the 
     * elements and changing elements with
     * partly real and partly imaginary elements to be real with
     * periodic real nodes instead of mirror image nodes.
     * 
     * @param rImageMap a map of four standard vectors of indices
     * 0. The left nodes which have been mirrored
     * 1. The image nodes relating to these left nodes (on right of mesh)
     * 2. The right nodes which have been mirrored
     * 3. The image nodes relating to these right nodes (on left of mesh)
     */
    void ReconstructCylindricalMesh(std::vector<std::vector<unsigned> >& rImageMap)
    {
        std::vector<unsigned> left_original = rImageMap[0];
        std::vector<unsigned> left_images = rImageMap[1];
        std::vector<unsigned> right_original = rImageMap[2];
        std::vector<unsigned> right_images = rImageMap[3];
        
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
                    if(IsThisIndexInList(this_node_index,left_images))
                    {
                        this_node_an_image = true;
                    }
                    if(IsThisIndexInList(this_node_index,right_images))
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
                 * Becuase of this we can just take one case (one image node)
                 * and delete the other element
                 */
                if (number_of_image_nodes==1)
                {
                    //std::cout << "Periodic element found \n" << std::flush;   
                    for (unsigned i=0 ; i<3 ; i++)
                    {
                        unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                        for (unsigned j=0 ; j<left_images.size() ; j++)
                        {
                            if(this_node_index==left_images[j])
                            {
                                p_element->ReplaceNode(mNodes[left_images[j]],mNodes[left_original[j]]);
                                //std::cout << "Node " << left_images[j] << " swapped for node " << left_original[j] << "\n" << std::flush;
                            }
                        }
                        for (unsigned j=0 ; j<right_images.size() ; j++)
                        {
                            if(this_node_index==right_images[j])
                            {
                                p_element->ReplaceNode(mNodes[right_images[j]],mNodes[right_original[j]]);
                                //std::cout << "Node " << right_images[j] << " swapped for node " << right_original[j] << "\n" << std::flush;
                            }
                        }
                    }
                }
            }
        }// end of loop over elements
        
        // Figure out which boundary elements have real nodes and image nodes in them
        // and replace image nodes with corresponding real ones.
        for (unsigned elem_index = 0; elem_index<GetNumAllBoundaryElements(); elem_index++)
        {
            BoundaryElement<1,2>* p_element = GetBoundaryElement(elem_index);
            if (!p_element->IsDeleted())
            {
                //std::cout << "Boundary Element " << elem_index << " connects nodes : ";
                unsigned number_of_image_nodes = 0;
                for (unsigned i=0 ; i<2 ; i++)
                {
                    unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                    //std::cout << this_node_index << "\t";
                    bool this_node_an_image = false;
                    if(IsThisIndexInList(this_node_index,left_images))
                    {
                        this_node_an_image = true;
                    }
                    if(IsThisIndexInList(this_node_index,right_images))
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
                    p_element->MarkAsDeleted();
                    mDeletedBoundaryElementIndices.push_back(p_element->GetIndex());
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
                        unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                        for (unsigned j=0 ; j<left_images.size() ; j++)
                        {
                            if(this_node_index==left_images[j])
                            {
                                //std::cout << "PERIODIC \n" << std::flush;  
                                p_element->ReplaceNode(mNodes[left_images[j]],mNodes[left_original[j]]);
                                //std::cout << "Node " << left_images[j] << " swapped for node " << left_original[j] << "\n" << std::flush;
                            }
                        }
                        for (unsigned j=0 ; j<right_images.size() ; j++)
                        {
                            if(this_node_index==right_images[j])
                            {
                                //std::cout << "IMAGE\n" << std::flush;
                                p_element->MarkAsDeleted();
                                mDeletedBoundaryElementIndices.push_back(p_element->GetIndex());
                            }
                        }
                    }
                }
            }
            
        }
        
        // Delete all image nodes
        for (unsigned i=0 ; i<left_images.size() ; i++)
        {
            mNodes[left_images[i]]->MarkAsDeleted();
            mDeletedNodeIndices.push_back(left_images[i]);
        }
        for (unsigned i=0 ; i<right_images.size() ; i++)
        {
            mNodes[right_images[i]]->MarkAsDeleted();
            mDeletedNodeIndices.push_back(right_images[i]);
        }
                
        // ReIndex the mesh
        ReIndex();
    }

    /**
     * This OVERWRITTEN method evaluates the (surface) distance between two points in a 2D Cylindrical geometry.
     * 
     * locations should lie between [mXleft, mXRight) 
     * 
     * @param rLocation1 the x and y co-ordinates of point 1
     * @param rLocation2 the x and y co-ordinates of point 2
     * 
     * @return the vector from location1 to location2
     */
    c_vector<double, 2> GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
    {
        assert(mXRight>mXLeft);
        assert(mXLeft<=rLocation1[0]);  // 1st point is not in cylinder
        assert(mXLeft<=rLocation2[0]);  // 2nd point is not in cylinder
        assert(mXRight>=rLocation1[0]);  // 1st point is not in cylinder
        assert(mXRight>=rLocation2[0]);  // 2nd point is not in cylinder
        
        c_vector<double, 2> vector;
        
        double two_pi_r = mXRight - mXLeft;
        
        double x1 = rLocation1[0];
        double x2 = rLocation2[0];
        double y1 = rLocation1[1];
        double y2 = rLocation2[1];
            
        double x_dist = x2 - x1;    // can be -ve
        double y_dist = y2 - y1;    // can be -ve
        
        // handle the cylindrical condition here
        // if the points are more than halfway around the cylinder apart
        // measure the other way.
        if ( x_dist > (two_pi_r / 2.0) )
        {
            x_dist = x_dist - two_pi_r;
        }
        if ( x_dist < -(two_pi_r / 2.0))
        {
            x_dist = x_dist + two_pi_r;  
        }
        
        vector[0] = x_dist;
        vector[1] = y_dist;
#define COVERAGE_IGNORE
        return vector;
#undef COVERAGE_IGNORE
    }
    
    /**
     * OVERWRITTEN function to set the location of a node.
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
    void SetNode(unsigned index, Point<2> point, bool concreteMove)
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
        if (point.rGetLocation()[0] >= mXRight)
        {   // move point to the left
            point.SetCoordinate(0u, point.rGetLocation()[0]-mWidth);
            //std::cout << "Moving point to the left\n" << std::flush;
        }
        if (point.rGetLocation()[0] < mXLeft)
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
    bool IsThisIndexInList(const unsigned& rNodeIndex, const std::vector<unsigned>& rListOfNodes)
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
    
    void TestTopAndBottomRowAlignment()
    {
        // Check that the top and bottom rows have the same y-co-ordinate 
        // or things will start to go wrong with boundary elements.
        double y_location = 0.0;
        for (unsigned i=0 ; i<mTopBoundary.size() ; i++)
        {
            y_location = mNodes[mTopBoundary[i]]->rGetLocation()[1];
            if (fabs(y_location - mTop)>1e-3)
            {
                //std::cout << "y = " << y_location << ", mTop = " << mTop << "\n" << std::flush;
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
     * OVERWRITTEN FUNCTION
     * @param rDimension must be 1 (x) or 2 (y)
     * @return width the CryptWidth or current height 
     */
    double GetWidth(const unsigned& rDimension)
    {
        double width=0.0;
        assert(rDimension==1 || rDimension==2);
        if (rDimension==1)
        {
            width = mWidth;   
        }
        if (rDimension==2)
        {
            width = mTop-mBottom;
        }
        return width;   
    }
    
    /**
     * OVERWRITTEN FUNCTION to ensure new node is introduced at a point on 
     * the cylinder and not slightly off it.
     *     
     * @param pNewNode pointer to a new node
     * @param map A node map of original mesh size
     * 
     * @return the new node index
     */
    unsigned AddNodeAndReMesh(Node<2> *pNewNode, NodeMap &map)
    {
        // Add the node to the mesh
        unsigned node_index = ConformingTetrahedralMesh<2,2>::AddNode(pNewNode);
        
        // If necessary move it to be back on the cylinder
        Point<2> new_node_point = pNewNode->GetPoint();
        SetNode(node_index, new_node_point, false);
        
        // increase the size of the node map to match the new mesh.
        map.Reserve(GetNumNodes());
        // Perform CYLINDRICAL ReMesh
        ReMesh(map);
        return node_index;
    }
};

#endif // _CYLINDRICAL2DMESH_CPP_


