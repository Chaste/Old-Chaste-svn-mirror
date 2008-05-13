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
#ifndef TESTCYLINDRICAL2DMESH_HPP_
#define TESTCYLINDRICAL2DMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <algorithm>

//#include "Cylindrical2dMesh.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TrianglesMeshWriter.hpp"

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
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true, crypt_width/cells_across);
        
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
        for (unsigned i=0; i<p_mesh->mLeftOriginals.size(); i++)
        {
            c_vector<double, 2> original_location = p_mesh->GetNode(p_mesh->mLeftOriginals[i])->rGetLocation();
            c_vector<double, 2> image_location = p_mesh->GetNode(p_mesh->mLeftImages[i])->rGetLocation();
            
            TS_ASSERT_DELTA(original_location[0]+crypt_width,image_location[0],1e-7);
            TS_ASSERT_DELTA(original_location[1],image_location[1],1e-7);
        }
        
        for (unsigned i=0; i<p_mesh->mRightOriginals.size(); i++)
        {
            c_vector<double, 2> original_location = p_mesh->GetNode(p_mesh->mRightOriginals[i])->rGetLocation();
            c_vector<double, 2> image_location = p_mesh->GetNode(p_mesh->mRightImages[i])->rGetLocation();
            
            TS_ASSERT_DELTA(original_location[0]-crypt_width,image_location[0],1e-7);
            TS_ASSERT_DELTA(original_location[1],image_location[1],1e-7);
        }
        
        // Check that we've got the correct number.
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up*2);
        
        // Cheat to put something into mTopHaloNodes to make exception throw.
        p_mesh->mTopHaloNodes.push_back(234u);
        unsigned corresponding_node_index = 0u;

        corresponding_node_index = p_mesh->GetCorrespondingNodeIndex(0u);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[0]+crypt_width,p_mesh->GetNode(corresponding_node_index)->rGetLocation()[0],1e-9 );
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[1],p_mesh->GetNode(corresponding_node_index)->rGetLocation()[1],1e-9);
        
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
//            ChastePoint<2> old_point(location);
//            c_vector<double, 2> new_location = location;
//            new_location[1] = -0.866025;
//            ChastePoint<2> new_point(new_location);
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
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
        // Create a mirrored load of nodes for the normal remesher to work with.
        p_mesh->CreateMirrorNodes();
        
        // Call the normal re-mesh
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ConformingTetrahedralMesh<2,2>::ReMesh(map);
        
        //
        // Re-Index the vectors regarding left/right nodes with the node map.
        //
        for (unsigned i = 0; i<p_mesh->mLeftOriginals.size(); i++)
        {
             p_mesh->mLeftOriginals[i]=map.GetNewIndex(p_mesh->mLeftOriginals[i]);
             p_mesh->mLeftImages[i]=map.GetNewIndex(p_mesh->mLeftImages[i]);
        }
        for (unsigned i = 0; i<p_mesh->mRightOriginals.size(); i++)
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
                for (unsigned i=0; i<3; i++)
                {
                    unsigned this_node_index = p_element->GetNodeGlobalIndex(i);
                    
                    if(this_node_index==0)
                    {   
                        elements_for_node_0++;
                        for (unsigned j=0; j<3; j++)
                        {
                            checksum_for_node_0 += p_element->GetNodeGlobalIndex(j);
                        }
                    }
                    if(this_node_index==11)
                    {   
                        elements_for_node_11++;
                        for (unsigned j=0; j<3; j++)
                        {
                            checksum_for_node_11 += p_element->GetNodeGlobalIndex(j);
                        }
                    }
                    if(this_node_index==12)
                    {   
                        elements_for_node_12++;
                        for (unsigned j=0; j<3; j++)
                        {
                            checksum_for_node_12 += p_element->GetNodeGlobalIndex(j);
                        }
                    }
                    if(this_node_index==18)
                    {   
                        elements_for_node_18++;
                        for (unsigned j=0; j<3; j++)
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
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
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
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
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
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
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
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
        c_vector<double, 2> location1 = p_mesh->GetNode(1)->rGetLocation();
        c_vector<double, 2> location2 = p_mesh->GetNode(4)->rGetLocation();

        // test a normal distance calculation
        c_vector<double, 2> vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 0.5, 1e-7);
        TS_ASSERT_DELTA(vector[1], sqrt(3.0)/2.0, 1e-4);
        // and the opposite vector
        vector = p_mesh->GetVectorFromAtoB(location2, location1);
        TS_ASSERT_DELTA(vector[0], -0.5, 1e-7);
        TS_ASSERT_DELTA(vector[1], -sqrt(3.0)/2.0, 1e-4);
        
        // test a periodic calculation
        location1[0] = 0.5;
        location1[1] = 3.0;
        location2[0] = 2.5;
        location2[1] = 4.0;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], -1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], +1.0, 1e-7);
        
        // test a periodic calculation where points need to be swapped
        location1[0] = 2.5;
        location1[1] = 4.0;
        location2[0] = 0.5;
        location2[1] = 3.0;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], +1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], -1.0, 1e-7);
        
        // we also want GetVectorFromAtoB to work when the x coord of A and B is not
        // between 0 and crypt width, by first normalizing the x coord (by taking modulus by crypt width)
        
        location1[0] = -0.5;
        location1[1] = 0;
        location2[0] = 2.5;
        location2[1] = 1;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 0.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], +1.0, 1e-7);
        
        // some tests where the location[0] is not between -0.25*crypt_width and 1.25*crypt_width
        
        location1[0] = -2.5;
        location1[1] = 0;
        location2[0] = 4.5;
        location2[1] = 1;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], +1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], +1.0, 1e-7);
        
    }
    
    void TestSetNodeLocationForCylindricalMesh() throw (Exception)
    {        
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        double crypt_width = 3.0;
        unsigned thickness_of_ghost_layer = 2;
        
        HoneycombMeshGenerator generator(cells_across,cells_up,thickness_of_ghost_layer,true,crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
        /*
         * New test to check that the top and bottom rows move together
         * bottom row = 0,1,2
         * top row = 18,19,20
         */
//        c_vector<double, 2> new_location = p_mesh->GetNode(0)->rGetLocation();
//        new_location[1] = -1.760;
//        ChastePoint<2> boundary_point(new_location);
//        // We just move one of the bottom boundary nodes and then...
//        p_mesh->SetNode(0u, boundary_point,false);
//        // check that all the nodes on this boundary have moved down
//        for (unsigned i=0; i<3; i++)
//        {
//            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetLocation()[1],-1.76000,1e-6);   
//        }
//        
//        // Same for one of the top boundary nodes
//        new_location = p_mesh->GetNode(19)->rGetLocation();
//        new_location[1] = 4.0;
//        ChastePoint<2> boundary_point2(new_location);
//        p_mesh->SetNode(19u, boundary_point2,false);
//        // check that all the nodes on this boundary have moved up
//        for (unsigned i=18; i<21; i++)
//        {
//            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetLocation()[1],4.0,1e-6);   
//        }
        
        // move one of the nodes to near the periodic boundary
        c_vector<double, 2> new_point_location;
        new_point_location[0] = 2.999999999;
        new_point_location[1] = -0.866025;
        ChastePoint<2> new_point(new_point_location);
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
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
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
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
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
        
        unsigned num_original_halo_nodes = p_mesh->GetNumNodes() - original_num_nodes;
        
        p_mesh->CreateMirrorNodes();
        
//        TrianglesMeshWriter<2,2> writer3("","HaloNodes.2");
//        writer3.WriteFilesUsingMesh(*p_mesh);
        
        double new_mesh_height = p_mesh->GetWidth(1);
        unsigned new_num_nodes = p_mesh->GetNumNodes();
        
        // Halo of nodes is added 0.5 above and below the original mesh.
        TS_ASSERT_DELTA(original_mesh_height, new_mesh_height, 1.0 + 1e-5);
        TS_ASSERT_EQUALS(new_num_nodes, original_num_nodes*2+2*num_original_halo_nodes);
        
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ConformingTetrahedralMesh<2,2>::ReMesh(map);   // recreates the boundary elements
        
//        TrianglesMeshWriter<2,2> writer4("","HaloNodes.3");
//        writer4.WriteFilesUsingMesh(*p_mesh);
        
        /*
         * TEST HALO NODE REMOVAL 
         */
        p_mesh->GenerateVectorsOfElementsStraddlingPeriodicBoundaries();
        p_mesh->CorrectNonPeriodicMesh();
        p_mesh->ReconstructCylindricalMesh();
        
//        TrianglesMeshWriter<2,2> writer5("","HaloNodes.4");
//        writer5.WriteFilesUsingMesh(*p_mesh);
        
        p_mesh->DeleteHaloNodes();
        
//        TrianglesMeshWriter<2,2> writer6("","HaloNodes.5");
//        writer6.WriteFilesUsingMesh(*p_mesh);
        
        TS_ASSERT_DELTA(original_mesh_height, p_mesh->GetWidth(1), 1.1e-6);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), original_num_nodes);
    }

    void TestHaloNodeReMesh() throw (Exception)
    {
        // This test checks that a Halo node remesh can handle a mesh of uneven height.
        
        unsigned cells_across = 5;
        unsigned cells_up = 3;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
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
        std::string dirname = "archive";
        OutputFileHandler handler(dirname, false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cylindrical_mesh_base.arch";
        
        std::string mesh_filename = "cylindrical_mesh";
        std::string mesh_pathname = handler.GetOutputDirectoryFullPath() + mesh_filename;
        
        // Set up a mesh
        unsigned cells_across = 5;
        unsigned cells_up = 3;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;
    
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true, crypt_width/cells_across);
        ConformingTetrahedralMesh<2,2> * const p_mesh = generator.GetCylindricalMesh();
        // You need the const above to stop a BOOST_STATIC_ASSERTION failure.
        // This is because the serialization library only allows you to save tracked
        // objects while the compiler considers them const, to prevent the objects changing
        // during the save, and so object tracking leading to wrong results.
        // (e.g. A is saved once via pointer, then changed, then saved again.  The second
        //  save notes that A was saved before, so doesn't write its data again, and the
        //  change is lost.)

        { // Serialize the mesh
            double width = p_mesh->GetWidth(0);
            TS_ASSERT_DELTA(width,crypt_width,1e-7);
            // Save the mesh data using mesh writers
            TrianglesMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
            mesh_writer.WriteFilesUsingMesh(*p_mesh);
            // Archive the mesh
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            // We have to serialize via a pointer here, or the derived class information is lost.
            output_arch << p_mesh;
        }
        
        { // De-serialize and compare
            ConformingTetrahedralMesh<2,2>* p_mesh2;
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> p_mesh2;
            // Re-initialise the mesh
            p_mesh2->Clear();
            TrianglesMeshReader<2,2> mesh_reader(mesh_pathname);
            p_mesh2->ConstructFromMeshReader(mesh_reader);
            
            // This is needed, as the honeycomb generator calls ReMesh prior to returning
            // the mesh, to make it 'properly cylindrical'.   However, the cylindrical ReMesh
            // creates lots of halo & mirror nodes prior to calling the base class ReMesh,
            // and doesn't clean up all the internal data structures when it removes them,
            // so some of the tests below fail, unless we also do a ReMesh after loading.
            // Even if a ReMesh isn't actually needed, it should be safe.
            NodeMap map(p_mesh2->GetNumNodes());
            p_mesh2->ReMesh(map);
            
            TS_ASSERT_DELTA(p_mesh2->GetWidth(0), crypt_width, 1e-7);
            
            // Compare the loaded mesh against the original
            TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(), p_mesh2->GetNumAllNodes());
            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), p_mesh2->GetNumNodes());
            TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryNodes(), p_mesh2->GetNumBoundaryNodes());
            
            for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
            {
                Node<2> *p_node = p_mesh->GetNode(i);
                Node<2> *p_node2 = p_mesh2->GetNode(i);
                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
                for (unsigned j=0; j<2; j++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-16);
                }
            }
            
            TS_ASSERT_EQUALS(p_mesh->GetNumElements(), p_mesh2->GetNumElements());
            TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), p_mesh2->GetNumAllElements());
            TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), p_mesh2->GetNumBoundaryElements());
            TS_ASSERT_EQUALS(p_mesh->GetNumAllBoundaryElements(), p_mesh2->GetNumAllBoundaryElements());
            ConformingTetrahedralMesh<2,2>::ElementIterator it=p_mesh->GetElementIteratorBegin();
            ConformingTetrahedralMesh<2,2>::ElementIterator it2=p_mesh2->GetElementIteratorBegin();
            for (;
                 it != p_mesh->GetElementIteratorEnd();
                 ++it, ++it2)
            {
                Element<2,2>* p_elt = *it;
                Element<2,2>* p_elt2 = *it2;
                TS_ASSERT_EQUALS(p_elt->GetNumNodes(), p_elt2->GetNumNodes());
                for (unsigned i=0; i<p_elt->GetNumNodes(); i++)
                {
                    TS_ASSERT_EQUALS(p_elt->GetNodeGlobalIndex(i), p_elt2->GetNodeGlobalIndex(i));
                }
            }
            
            // We now need to free the mesh, since there is no honeycomb generator to do so.
            delete p_mesh2;
        }
    }
    
    void TestConstructFromNodeList()
    {
        std::vector<Node<2> *> nodes;
        
        nodes.push_back(new Node<2>(0, true,  0.1,  -0.01));
        nodes.push_back(new Node<2>(1, true,  0.5,  0.0));
        nodes.push_back(new Node<2>(2, true,  0.9,  0.01));
        nodes.push_back(new Node<2>(3, true,  0.1,  0.99));
        nodes.push_back(new Node<2>(4, true,  0.5,  1.0));
        nodes.push_back(new Node<2>(5, true,  0.9,  1.01));
        
        const double width = 1.0;
        Cylindrical2dMesh mesh(width, nodes);
        
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6u);
        
        // find the element with node indices 2,3,5 (which stradles the periodic boundary)
        unsigned element_index;
        std::set<unsigned> target_element_node_indices;
        for (element_index=0; element_index<mesh.GetNumElements(); element_index++)
        {
            target_element_node_indices.clear();
            target_element_node_indices.insert(2);
            target_element_node_indices.insert(3);
            target_element_node_indices.insert(5);
            
            for (unsigned node_local_index = 0; node_local_index<=2; node_local_index++)
            {
                target_element_node_indices.erase(mesh.GetElement(element_index)->GetNodeGlobalIndex(node_local_index));
            }
            if (target_element_node_indices.empty())
            {
                break;
            }
        }
        TS_ASSERT(target_element_node_indices.empty());
        
        
        // calculate the circumsphere of the element
        c_vector<double, 3> circumsphere = mesh.GetElement(element_index)->CalculateCircumsphere();
        
        TS_ASSERT_DELTA(circumsphere[0], 0.9509, 1e-3);
        TS_ASSERT_DELTA(circumsphere[1], 0.5100, 1e-3);
        TS_ASSERT_DELTA(circumsphere[2], 0.2526, 1e-3);

        /* The reason that the circumsphere is calculated correctly for a periodic boundary
         * stradling element is somewhat obscure
         * The jacobian of the element is calculate when the element has a mirror node
         * The mirror node is then replaced with the node within the periodic mesh
         * The circumsphere is calculated based on the jacobian and the replaced node within the periodic mesh
         * 
         * uncommenting the following line of code causes an error:
         * Jacobian determinant is non-positive
         * 
         * mesh.GetElement(element_index)->RefreshJacobianDeterminant();
         */
    }
    
    void TestGenerateVectorsOfElementsStraddlingPeriodicBoundaries()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/bad_cylindrical_9_1");
        Cylindrical2dMesh mesh(9.1);
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // we now emulate the commands of the ReMesh function as far as it goes before Generating the lists.
        {
            mesh.CreateHaloNodes();
            mesh.CreateMirrorNodes();
            
            NodeMap big_map(mesh.GetNumAllNodes()); 
            mesh.ConformingTetrahedralMesh<2,2>::ReMesh(big_map);
        }
        
        mesh.GenerateVectorsOfElementsStraddlingPeriodicBoundaries();
    
        TS_ASSERT_EQUALS(mesh.mLeftPeriodicBoundaryElementIndices.size(), 43u);
        TS_ASSERT_EQUALS(mesh.mRightPeriodicBoundaryElementIndices.size(), 42u);
        
        // TEST the GetCorrespondingNodeIndex() method
        TS_ASSERT_EQUALS(mesh.GetCorrespondingNodeIndex(393), 84u);
        TS_ASSERT_EQUALS(mesh.GetCorrespondingNodeIndex(84), 393u);
        TS_ASSERT_EQUALS(mesh.GetCorrespondingNodeIndex(188), 329u);
        TS_ASSERT_EQUALS(mesh.GetCorrespondingNodeIndex(329), 188u);
        
        mesh.CorrectNonPeriodicMesh();

        mesh.DeleteHaloNodes();
    }
    
    void TestCorrectNonPeriodicMeshMapLeftToRight()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/bad_cylindrical_9_1");
        Cylindrical2dMesh mesh(9.1);
        mesh.ConstructFromMeshReader(mesh_reader);
        
        NodeMap map(0);
        mesh.ReMesh(map);
        assert(map.IsIdentityMap());
        
        for (unsigned node_index=0; node_index<mesh.GetNumAllNodes(); node_index++)
        {
            std::vector<unsigned> indices;
            
            
            //Get the forward star from each node that isn't at the top or bottom boundary
            Node<2> *p_node=mesh.GetNode(node_index);
            if (p_node->rGetLocation()[1] < -2.5)
            {
                continue;
            }
            if (p_node->rGetLocation()[1] > 13.8)
            {
                continue;
            }
                
            //Iterate over countaining elements to get the elements of the forward star
            for (Node<2>::ContainingElementIterator it = p_node->ContainingElementsBegin();
                it != p_node->ContainingElementsEnd();
                ++it)
            {
                Element <2, 2> *p_element = mesh.GetElement(*it);
                for (unsigned j=0; j<3; j++)
                {
                    unsigned index=p_element->GetNodeGlobalIndex(j);
                    if (index != node_index)
                    {
                        indices.push_back(index);
                    }
                }
                
                
            }
            //Each node in the forward star should appear exactly twice.  Sort and test.
            sort(indices.begin(), indices.end());
            for (unsigned i=0; i<indices.size();i++)
            {
                if (i%2 == 0)
                {
                   TS_ASSERT_EQUALS(indices[i],indices[i+1])
                }
                
            }
                 
        }
            
    }        
        
};

#endif /*TESTCYLINDRICAL2DMESH_HPP_*/


