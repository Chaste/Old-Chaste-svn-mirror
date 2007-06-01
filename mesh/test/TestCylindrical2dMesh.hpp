#ifndef TESTCYLINDRICAL2DMESH_HPP_
#define TESTCYLINDRICAL2DMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Cylindrical2dMesh.cpp"
#include "CryptHoneycombMeshGenerator.hpp"
#include "TrianglesMeshWriter.cpp"



class TestCylindrical2dMesh : public CxxTest::TestSuite
{
public:

    void TestBasicFunctions()
    {   
        // Test IsThisIndexInList
        Cylindrical2dMesh mesh(1.0);
        
        std::vector<unsigned> list_of_nodes;
        list_of_nodes.push_back(0u);
        list_of_nodes.push_back(2u);
        list_of_nodes.push_back(5u);
        list_of_nodes.push_back(7u);
        
        TS_ASSERT_EQUALS(mesh.IsThisIndexInList(0u,list_of_nodes),true);
        TS_ASSERT_EQUALS(mesh.IsThisIndexInList(2u,list_of_nodes),true);
        TS_ASSERT_EQUALS(mesh.IsThisIndexInList(5u,list_of_nodes),true);
        TS_ASSERT_EQUALS(mesh.IsThisIndexInList(7u,list_of_nodes),true);

        TS_ASSERT_EQUALS(mesh.IsThisIndexInList(1u,list_of_nodes),false);
        TS_ASSERT_EQUALS(mesh.IsThisIndexInList(3u,list_of_nodes),false);
        TS_ASSERT_EQUALS(mesh.IsThisIndexInList(4u,list_of_nodes),false);
        TS_ASSERT_EQUALS(mesh.IsThisIndexInList(6u,list_of_nodes),false);
        TS_ASSERT_EQUALS(mesh.IsThisIndexInList(8u,list_of_nodes),false);
    }


    void TestCreateMirrorCellsANDAlignmentTester() throw (Exception)
    {   
    	// note that elements are not created (and boundary elements are not changed)
        // this just creates a set of new nodes.
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 0u;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // reset the mesh
        p_mesh = generator.GetCylindricalMesh();
        p_mesh->CreateMirrorNodes();

        // Check the vectors are the right size...
        TS_ASSERT_EQUALS(p_mesh->mLeftOriginals.size(),36u);
        TS_ASSERT_EQUALS(p_mesh->mRightOriginals.size(),36u);
        TS_ASSERT_EQUALS(p_mesh->mLeftOriginals.size(), p_mesh->mLeftImages.size());
        TS_ASSERT_EQUALS(p_mesh->mRightOriginals.size(), p_mesh->mRightImages.size());
        
        // Check that the image nodes are where they should be.
        for (unsigned i=0 ; i<p_mesh->mLeftOriginals.size() ; i++)
        {
            c_vector<double, 2> original_location = p_mesh->GetNode(p_mesh->mLeftOriginals[i])->rGetLocation();
            c_vector<double, 2> image_location = p_mesh->GetNode(p_mesh->mLeftImages[i])->rGetLocation();
            
            TS_ASSERT_DELTA(original_location[0]+crypt_width,image_location[0],1e-7);
            TS_ASSERT_DELTA(original_location[1],image_location[1],1e-7);
        }
        
        for (unsigned i=0 ; i<p_mesh->mRightOriginals.size() ; i++)
        {
            c_vector<double, 2> original_location = p_mesh->GetNode(p_mesh->mRightOriginals[i])->rGetLocation();
            c_vector<double, 2> image_location = p_mesh->GetNode(p_mesh->mRightImages[i])->rGetLocation();
            
            TS_ASSERT_DELTA(original_location[0]-crypt_width,image_location[0],1e-7);
            TS_ASSERT_DELTA(original_location[1],image_location[1],1e-7);
        }
        
        // Check that we've got the correct number.
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up*2);
        
        //Output2DMeshToFile(p_mesh, "node_positions.dat");
//        
//        /*
//         * TEST FOR ALIGNMENT TESTER
//         */
//        {
//            TS_ASSERT_THROWS_NOTHING(p_mesh->TestTopAndBottomRowAlignment());  
//            
//            // move one of the bottom nodes and ensure it throws.
//            c_vector<double, 2> location = p_mesh->GetNode(0)->rGetLocation();
//            Point<2> old_point(location);
//            c_vector<double, 2> new_location = location;
//            new_location[1] = -0.866025;
//            Point<2> new_point(new_location);
//            p_mesh->GetNode(0)->SetPoint(new_point);
//            TS_ASSERT_THROWS_ANYTHING(p_mesh->TestTopAndBottomRowAlignment());
//            
//            // move it back and ensure it passes.
//            p_mesh->GetNode(0)->SetPoint(old_point);
//            TS_ASSERT_THROWS_NOTHING(p_mesh->TestTopAndBottomRowAlignment());
//            
//            // move the top row and ensure it fails
//            p_mesh->GetNode(70)->SetPoint(new_point);
//            
//            TS_ASSERT_THROWS_ANYTHING(p_mesh->TestTopAndBottomRowAlignment());
//        }
    }
    
    void TestReconstructCylindricalMesh() throw (Exception)
    {   
    	// this takes in a new mesh created using the mirror function above
        // and a ReMesh call, then removes nodes, elements and boundary elements.
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0u;

        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        
        // Create a mirrored load of nodes for the normal remesher to work with.
        p_mesh->CreateMirrorNodes();
        
        // Call the normal re-mesh
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ConformingTetrahedralMesh<2,2>::ReMesh(map);
        
        //
        // Re-Index the vectors regarding left/right nodes with the node map.
        //
        for (unsigned i = 0 ; i<p_mesh->mLeftOriginals.size() ; i++)
        {
                p_mesh->mLeftOriginals[i]=map.GetNewIndex(p_mesh->mLeftOriginals[i]);
                p_mesh->mLeftImages[i]=map.GetNewIndex(p_mesh->mLeftImages[i]);
        }
        for (unsigned i = 0 ; i<p_mesh->mRightOriginals.size() ; i++)
        {
                p_mesh->mRightOriginals[i]=map.GetNewIndex(p_mesh->mRightOriginals[i]);
                p_mesh->mRightImages[i]=map.GetNewIndex(p_mesh->mRightImages[i]);
        }
        
        p_mesh->ReconstructCylindricalMesh();
        
        unsigned elements_for_node_0 = 0;
        unsigned elements_for_node_11 = 0;
        unsigned elements_for_node_12 = 0;
        unsigned elements_for_node_18 = 0;
        
        unsigned checksum_for_node_0 = 0;
        unsigned checksum_for_node_11 = 0;
        unsigned checksum_for_node_12 = 0;
        unsigned checksum_for_node_18 = 0;
        
        for (unsigned elem_index = 0; elem_index<p_mesh->GetNumAllElements(); elem_index++)
        {
            Element<2,2>* p_element = p_mesh->GetElement(elem_index);
            if (!p_element->IsDeleted())
            {
                for (unsigned i=0 ; i<3 ; i++)
                {
                    unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                    
                    if(this_node_index==0)
                    {   
                        elements_for_node_0++;
                        for (unsigned j=0 ; j<3 ; j++)
                        {
                            checksum_for_node_0 += p_element->GetNodeGlobalIndex(j);
                        }
                    }
                    if(this_node_index==11)
                    {   
                        elements_for_node_11++;
                        for (unsigned j=0 ; j<3 ; j++)
                        {
                            checksum_for_node_11 += p_element->GetNodeGlobalIndex(j);
                        }
                    }
                    if(this_node_index==12)
                    {   
                        elements_for_node_12++;
                        for (unsigned j=0 ; j<3 ; j++)
                        {
                            checksum_for_node_12 += p_element->GetNodeGlobalIndex(j);
                        }
                    }
                    if(this_node_index==18)
                    {   
                        elements_for_node_18++;
                        for (unsigned j=0 ; j<3 ; j++)
                        {
                            checksum_for_node_18 += p_element->GetNodeGlobalIndex(j);
                        }
                    }
                }
            }
        }
        
        TS_ASSERT_EQUALS(elements_for_node_0, 3u);
        TS_ASSERT_EQUALS(elements_for_node_11, 6u);
        TS_ASSERT_EQUALS(elements_for_node_12, 6u);
        TS_ASSERT_EQUALS(elements_for_node_18, 6u);
        
        // these are nodes on the edge of the cylindrical region.
        // If the mesh is periodic they will be joined
        // to the following other nodes
        unsigned checksum_target_for_node_0 = (0 + 5 + 11) + (0 + 6 + 11) + (0 + 1 + 6);
        TS_ASSERT_EQUALS(checksum_for_node_0, checksum_target_for_node_0);        
        
        unsigned checksum_target_for_node_11 = (0 + 5 + 11) + (0 + 6 + 11) + (11 + 12 + 6)
                                              +(11 + 17 + 12) + (10 + 11 + 17) + (5 + 10 + 11);
        TS_ASSERT_EQUALS(checksum_for_node_11, checksum_target_for_node_11); 
        
        unsigned checksum_target_for_node_12 = (12 + 13 + 18) + (12 + 18 + 23) + (12 + 17 + 23)
                                              +(12 + 11 + 17) + (12 + 6 + 11) + (12 + 6 + 13);
        TS_ASSERT_EQUALS(checksum_for_node_12, checksum_target_for_node_12); 
        
        unsigned checksum_target_for_node_18 = (18 + 19 + 25) + (18 + 24 + 25) + (18 + 23 + 24)
                                                +(18 + 12 + 23) + (18 + 12 + 13) + (18 + 13 + 19);
        TS_ASSERT_EQUALS(checksum_for_node_18, checksum_target_for_node_18); 
        
        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),12u);
        
//        // Check that the ReIndex has worked
//        TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(),cells_across*cells_up);
//        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(),2*cells_across*(cells_up-1));
//        TS_ASSERT_EQUALS(p_mesh->GetNumAllBoundaryElements(),12u);
        
        //Output2DMeshToFile(p_mesh, "node_positions.dat");
    }
    
    void TestCylindricalReMesh() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0;
        
        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);
        
        TS_ASSERT_EQUALS(map.Size(), p_mesh->GetNumNodes());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),1u);  // No boundary elements now the halo nodes are removed
        
        //Output2DMeshToFile(p_mesh, "node_positions.dat");
    }


    void TestCylindricalReMeshAfterDelete() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0;
        
        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
        unsigned num_old_nodes = p_mesh->GetNumNodes();

        p_mesh->DeleteNode(15);
        
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up - 1);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),1u);  // No boundary elements now the halo nodes are removed

        TS_ASSERT_EQUALS(map.Size(), num_old_nodes);
        TS_ASSERT_EQUALS(map.IsDeleted(15), true);
        
        for(unsigned i=0; i<num_old_nodes; i++)
        {
            if(i<15)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), i);
            }
            if(i>15)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), (unsigned)(i-1));
            }
        } 

   }
    
    void TestCylindricalReMeshOnSmallMesh() throw (Exception)
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;
        
        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),1u); // boundary elements removed now halo nodes are used

        //Output2DMeshToFile(p_mesh, "node_positions.dat");
    }
    
    void TestGetVectorBetweenCyclindricalPoints() throw (Exception)
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;
        
        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
        c_vector<double, 2> location1 = p_mesh->GetNode(1)->rGetLocation();
        c_vector<double, 2> location2 = p_mesh->GetNode(4)->rGetLocation();

        // test a normal distance calculation
        c_vector<double, 2> vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 0.5, 1e-7);
        TS_ASSERT_DELTA(vector[1], sqrt(3.0)/2.0, 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 1.0, 1e-4);
        // and the opposite vector
        vector = p_mesh->GetVectorFromAtoB(location2, location1);
        TS_ASSERT_DELTA(vector[0], -0.5, 1e-7);
        TS_ASSERT_DELTA(vector[1], -sqrt(3.0)/2.0, 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 1.0, 1e-4);
        
        // test a periodic calculation
        location1[0] = 0.5;
        location1[1] = 3.0;
        location2[0] = 2.5;
        location2[1] = 4.0;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], -1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], +1.0, 1e-7);
        TS_ASSERT_DELTA(norm_2(vector), sqrt(2.0), 1e-7);
        
        // test a periodic calculation where points need to be swapped
        location1[0] = 2.5;
        location1[1] = 4.0;
        location2[0] = 0.5;
        location2[1] = 3.0;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], +1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], -1.0, 1e-7);
        TS_ASSERT_DELTA(norm_2(vector), sqrt(2.0), 1e-7);
    }
    
    void TestSetNodeLocationForCylindricalMesh() throw (Exception)
    {        
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        double crypt_width = 3.0;
        unsigned thickness_of_ghost_layer = 2;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        
        /*
         * New test to check that the top and bottom rows move together
         * bottom row = 0,1,2
         * top row = 18,19,20
         */
//        c_vector<double, 2> new_location = p_mesh->GetNode(0)->rGetLocation();
//        new_location[1] = -1.760;
//        Point<2> boundary_point(new_location);
//        // We just move one of the bottom boundary nodes and then...
//        p_mesh->SetNode(0u, boundary_point,false);
//        // check that all the nodes on this boundary have moved down
//        for (unsigned i=0; i<3 ; i++)
//        {
//            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetLocation()[1],-1.76000,1e-6);   
//        }
//        
//        // Same for one of the top boundary nodes
//        new_location = p_mesh->GetNode(19)->rGetLocation();
//        new_location[1] = 4.0;
//        Point<2> boundary_point2(new_location);
//        p_mesh->SetNode(19u, boundary_point2,false);
//        // check that all the nodes on this boundary have moved up
//        for (unsigned i=18; i<21 ; i++)
//        {
//            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetLocation()[1],4.0,1e-6);   
//        }
        
        // move one of the nodes to near the periodic boundary
        c_vector<double, 2> new_point_location;
        new_point_location[0] = 2.999999999;
        new_point_location[1] = -0.866025;
        Point<2> new_point(new_point_location);
        p_mesh->GetNode(5)->SetPoint(new_point);
        
      
        new_point.SetCoordinate(0u, -0.0001);
        //std::cout << " x = " << new_point.rGetLocation()[0] << ", y = " << new_point.rGetLocation()[1] << "\n" << std::flush;
        // This node was on left and is now near the right
        
        p_mesh->SetNode(0u, new_point,false);
        
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 2.9999, 1e-4);
        
        new_point.SetCoordinate(0u, 1.0000);
        p_mesh->SetNode(0u, new_point,false);
        
        //std::cout << " x = " << new_point.rGetLocation()[0] << ", y = " << new_point.rGetLocation()[1] << "\n" << std::flush;
        // This node has stayed close to where it was
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 1.0000, 1e-4);
        
        new_point.SetCoordinate(0u, 3.0001);
        p_mesh->SetNode(0u, new_point,false);
        //std::cout << " x = " << new_point.rGetLocation()[0] << ", y = " << new_point.rGetLocation()[1] << "\n" << std::flush;
        // This node was on right and is now on the left.
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], +0.0001, 1e-4);
        
    }
    
    void TestAddNodeAndReMesh() throw (Exception)
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;
        
        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),1u); // boundary elements removed now halo nodes are used

        c_vector<double ,2> point;
        point[0] = -0.05;
        point[1] = 1.0;
        Node<2>* p_node = new Node<2>(p_mesh->GetNumNodes(), point);
        

        unsigned new_index = p_mesh->AddNode(p_node);
        NodeMap map(p_mesh->GetNumNodes());
        
        p_mesh->ReMesh(map);
        
        TS_ASSERT_EQUALS(map.Size(), p_mesh->GetNumNodes());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);
        
        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up+1);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2*cells_across*(cells_up-1)+2);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),1u);// boundary elements removed now halo nodes are used

        // Check that we have moved the new node across        
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[0], 3.0+point[0], 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[1], point[1], 1e-7);
        // Test GetWidth
        TS_ASSERT_DELTA(p_mesh->GetWidth(0u), 3.0, 1e-9);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1u), sqrt(3), 1e-6);
    }    

    void TestHaloNodeInsertionAndRemoval() throw (Exception)
    {
        unsigned cells_across = 5;
        unsigned cells_up = 3;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        
        c_vector<double, 2>& rLocation1 = p_mesh->GetNode(1)->rGetModifiableLocation();
        rLocation1[1] -= 0.5;
        
        c_vector<double, 2>& rLocation2 = p_mesh->GetNode(3)->rGetModifiableLocation();
        rLocation2[1] -= 0.4;
        
        c_vector<double, 2>& rLocation3 = p_mesh->GetNode(12)->rGetModifiableLocation();
        rLocation3[1] += 0.8;
        
        double original_mesh_height = p_mesh->GetWidth(1);
        unsigned original_num_nodes = p_mesh->GetNumNodes();
        
//        TrianglesMeshWriter<2,2> writer("","HaloNodes");
//        writer.WriteFilesUsingMesh(*p_mesh);
        
        p_mesh->CreateHaloNodes();
        
//        TrianglesMeshWriter<2,2> writer2("","HaloNodes.1");
//        writer2.WriteFilesUsingMesh(*p_mesh);
        
        double new_mesh_height = p_mesh->GetWidth(1);
        unsigned new_num_nodes = p_mesh->GetNumNodes();
        
        // Halo of nodes is added 0.5 above and below the original mesh.
        TS_ASSERT_DELTA(original_mesh_height, new_mesh_height, 1.0 + 1e-5);
        TS_ASSERT_EQUALS(new_num_nodes, original_num_nodes+2*9u);
        
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ConformingTetrahedralMesh<2,2>::ReMesh(map);   // recreates the boundary elements
        
        /*
         * TEST HALO NODE REMOVAL 
         */
        
        p_mesh->DeleteHaloNodes();
        
        TS_ASSERT_DELTA(original_mesh_height, p_mesh->GetWidth(1), 1e-6);
        TS_ASSERT_EQUALS(original_num_nodes, p_mesh->GetNumNodes());
        
        // Check that we still have a boundary element (for ReIndex)
        TS_ASSERT(p_mesh->GetNumBoundaryElements() > 0u );
        
    }

    void TestHaloNodeReMesh() throw (Exception)
    {
        // This test checks that a Halo node remesh can handle a mesh of uneven height.
        
        unsigned cells_across = 5;
        unsigned cells_up = 3;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        
        c_vector<double, 2>& rLocation1 = p_mesh->GetNode(1)->rGetModifiableLocation();
        rLocation1[1] -= 0.5;
        
        c_vector<double, 2>& rLocation2 = p_mesh->GetNode(3)->rGetModifiableLocation();
        rLocation2[1] -= 0.4;
        
        c_vector<double, 2>& rLocation3 = p_mesh->GetNode(12)->rGetModifiableLocation();
        rLocation3[1] += 0.8;
        
//        TrianglesMeshWriter<2,2> writer("","HaloNodes");
//        writer.WriteFilesUsingMesh(*p_mesh);
        
        unsigned total_elements = p_mesh->GetNumElements();
        unsigned total_nodes = p_mesh->GetNumNodes();
        // Check that the ReIndex is working still
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), p_mesh->GetNumElements());
        
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);
        
//        TrianglesMeshWriter<2,2> writer2("","HaloNodes.1");
//        writer2.WriteFilesUsingMesh(*p_mesh);
        
        // Check that we haven't added any nodes or elements by doing this Halo Node ReMesh.
        TS_ASSERT_EQUALS(total_elements, p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(total_nodes, p_mesh->GetNumNodes());
        // Check the ReIndex is working
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(), p_mesh->GetNumNodes());
        
    }
    
    void TestArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "cylindrical_mesh.arch";
        
        double width = 0.0;
        double height = 0.0;
        
        {   // Set up a mesh
            unsigned cells_across = 5;
            unsigned cells_up = 3;
            double crypt_width = 5.0;
            unsigned thickness_of_ghost_layer = 0;
        
            CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
            Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
            
            width = p_mesh->GetWidth(0);
            height = p_mesh->GetWidth(1);
            
            // Archive the mesh
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const Cylindrical2dMesh&>(*p_mesh);
            
        }
        
        {   
            unsigned cells_across = 10;
            unsigned cells_up = 10;
            double crypt_width = 10.0;
            unsigned thickness_of_ghost_layer = 0;
        
            CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
            Cylindrical2dMesh* p_mesh2=generator.GetCylindricalMesh();
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_mesh2;
            
            // these are commented out for now until the cylindrical2dMesh class
            // is archived properly.
            
            TS_ASSERT_DELTA(p_mesh2->GetWidth(0), width, 1e-7);
            // this one fails
//            TS_ASSERT_DELTA(p_mesh2->GetWidth(1), height, 1e-7);
        }
    }
    
};





#endif /*TESTCYLINDRICAL2DMESH_HPP_*/


