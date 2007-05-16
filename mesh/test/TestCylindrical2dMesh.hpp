#ifndef TESTCYLINDRICAL2DMESH_HPP_
#define TESTCYLINDRICAL2DMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "Cylindrical2dMesh.cpp"
#include "CryptHoneycombMeshGenerator.hpp"
#include "TrianglesMeshWriter.cpp"


class TestCylindrical2dMesh : public CxxTest::TestSuite
{
public:

    void TestBasicFunctions()
    {
        // Test default constructor cannot be used.
        TS_ASSERT_THROWS_ANYTHING(Cylindrical2dMesh mesh);
        
        std::vector <unsigned> empty;
        
        // Test IsThisIndexInList
        Cylindrical2dMesh mesh(0.0, 1.0, 2.0, -1.0, empty,empty);
        
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
        
        // Test GetWidth
        TS_ASSERT_DELTA(mesh.GetWidth(1u), 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetWidth(2u), 3.0, 1e-9);
        
    }

    void TestCreateMirrorCellsANDAlignmentTester() throw (Exception)
    {   // note that elements are not created (and boundary elements are not changed)
        // this just creates a set of new nodes.
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 0u;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
        std::vector<std::vector<unsigned> > image_map;

        // reset the mesh
        p_mesh = generator.GetCylindricalMesh();
        image_map = p_mesh->CreateMirrorNodes();
        
        //Output2DMeshToFile(p_mesh, "node_positions.dat");
        
        std::vector<unsigned> left_original = image_map[0];
        std::vector<unsigned> left_images = image_map[1];
        std::vector<unsigned> right_original = image_map[2];
        std::vector<unsigned> right_images = image_map[3];
        
        //std::cout << "Left size = " << left_original.size()<< ", Right size = " << right_original.size()<< "\n" << std::flush;
        // Check the vectors are the right size...
        TS_ASSERT_EQUALS(left_original.size(),36u);
        TS_ASSERT_EQUALS(right_original.size(),36u);
        TS_ASSERT_EQUALS(left_original.size(), left_images.size());
        TS_ASSERT_EQUALS(right_original.size(), right_images.size());
        
        // Check that the image nodes are where they should be.
        for (unsigned i=0 ; i<left_original.size() ; i++)
        {
            c_vector<double, 2> original_location = p_mesh->GetNode(left_original[i])->rGetLocation();
            c_vector<double, 2> image_location = p_mesh->GetNode(left_images[i])->rGetLocation();
            
            TS_ASSERT_DELTA(original_location[0]+crypt_width,image_location[0],1e-7);
            TS_ASSERT_DELTA(original_location[1],image_location[1],1e-7);
        }
        
        for (unsigned i=0 ; i<right_original.size() ; i++)
        {
            c_vector<double, 2> original_location = p_mesh->GetNode(right_original[i])->rGetLocation();
            c_vector<double, 2> image_location = p_mesh->GetNode(right_images[i])->rGetLocation();
            
            TS_ASSERT_DELTA(original_location[0]-crypt_width,image_location[0],1e-7);
            TS_ASSERT_DELTA(original_location[1],image_location[1],1e-7);
        }
        
        // Check that we've got the correct number.
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up*2);
        
        //Output2DMeshToFile(p_mesh, "node_positions.dat");
        
        /*
         * TEST FOR ALIGNMENT TESTER
         */
        {
            TS_ASSERT_THROWS_NOTHING(p_mesh->TestTopAndBottomRowAlignment());  
            
            // move one of the bottom nodes and ensure it throws.
            c_vector<double, 2> location = p_mesh->GetNode(0)->rGetLocation();
            Point<2> old_point(location);
            c_vector<double, 2> new_location = location;
            new_location[1] = -0.866025;
            Point<2> new_point(new_location);
            p_mesh->GetNode(0)->SetPoint(new_point);
            TS_ASSERT_THROWS_ANYTHING(p_mesh->TestTopAndBottomRowAlignment());
            
            // move it back and ensure it passes.
            p_mesh->GetNode(0)->SetPoint(old_point);
            TS_ASSERT_THROWS_NOTHING(p_mesh->TestTopAndBottomRowAlignment());
            
            // move the top row and ensure it fails
            p_mesh->GetNode(70)->SetPoint(new_point);
            
            TS_ASSERT_THROWS_ANYTHING(p_mesh->TestTopAndBottomRowAlignment());
        }
    }
    
    void TestReconstructCylindricalMesh() throw (Exception)
    {   // this takes in a new mesh created using the mirror function above
        // and a ReMesh call, then removes nodes, elements and boundary elements.
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0u;

        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        
        // Create a mirrored load of nodes for the normal remesher to work with.
        std::vector<std::vector<unsigned> > image_map = p_mesh->CreateMirrorNodes();
    
        std::vector<unsigned> left_original = image_map[0];
        std::vector<unsigned> left_images = image_map[1];
        std::vector<unsigned> right_original = image_map[2];
        std::vector<unsigned> right_images = image_map[3];
        
        // Call the normal re-mesh
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ConformingTetrahedralMesh<2,2>::ReMesh(map);
        
        //
        // Re-Index the image_map according to the node_map.
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
        
        image_map[0] = left_original;
        image_map[1] = left_images;
        image_map[2] = right_original;
        image_map[3] = right_images;
        
        p_mesh->ReconstructCylindricalMesh(image_map);
        
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
        
        // Check that the ReIndex has worked
        TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(),cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(),2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumAllBoundaryElements(),12u);
        
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

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),2*cells_across);
        
        //Output2DMeshToFile(p_mesh, "node_positions.dat");
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
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),2*cells_across);

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
        c_vector<double, 2> new_location = p_mesh->GetNode(0)->rGetLocation();
        new_location[1] = -1.760;
        Point<2> boundary_point(new_location);
        // We just move one of the bottom boundary nodes and then...
        p_mesh->SetNode(0u, boundary_point,false);
        // check that all the nodes on this boundary have moved down
        for (unsigned i=0; i<3 ; i++)
        {
            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetLocation()[1],-1.76000,1e-6);   
        }
        
        // Same for one of the top boundary nodes
        new_location = p_mesh->GetNode(19)->rGetLocation();
        new_location[1] = 4.0;
        Point<2> boundary_point2(new_location);
        p_mesh->SetNode(19u, boundary_point2,false);
        // check that all the nodes on this boundary have moved up
        for (unsigned i=18; i<21 ; i++)
        {
            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetLocation()[1],4.0,1e-6);   
        }
        
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
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),2*cells_across);

        c_vector<double ,2> point;
        point[0] = -0.05;
        point[1] = 1.0;
        Node<2>* p_node = new Node<2>(p_mesh->GetNumNodes(), point);
        
        NodeMap map(p_mesh->GetNumNodes());
        unsigned new_index = p_mesh->AddNode(p_node);
        p_mesh->ReMesh(map);
        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up+1);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2*cells_across*(cells_up-1)+2);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),2*cells_across);

        // Check that we have moved the new node across        
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[0], 3.0+point[0], 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[1], point[1], 1e-7);

        
    }    

};





#endif /*TESTCYLINDRICAL2DMESH_HPP_*/


