#ifndef TESTREMESH_HPP_
#define TESTREMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "OutputFileHandler.hpp"
#include "CryptHoneycombMeshGenerator.hpp"
#include "CancerParameters.hpp"
#include "NodeMap.hpp"
#include <iostream>
#include <cmath>

#include <vector>

class TestRemesh : public CxxTest::TestSuite
{
private:
	void Output2DMeshToFile(ConformingTetrahedralMesh<2,2>* p_mesh, std::string fileName)
    {
    	OutputFileHandler handler("");
    	out_stream file=handler.OpenOutputFile(fileName);
    	
    	unsigned num_nodes=p_mesh->GetNumNodes();
    
		for (unsigned i=0; i<num_nodes; i++)
    	{
        	c_vector<double, 2> location = p_mesh->GetNode(i)->rGetLocation();
        	(*file) << location[0] << "\t" << location[1] << "\n" << std::flush;
    	}
    	
    	file->close();
    }
	
public:
    void TestOperationOfTriangle() throw (Exception)
    {
        OutputFileHandler handler("");
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/SquarePartDecimation");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        double area=mesh.CalculateMeshVolume();
        TS_ASSERT_DELTA(0.01, area, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),77U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(),141U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),11U);
        
        out_stream node_file=handler.OpenOutputFile("temp.node");
        (*node_file)<<mesh.GetNumNodes()<<"\t2\t0\t0\n";
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            const c_vector<double, 2> node_loc = mesh.GetNode(i)->rGetLocation();
            (*node_file)<<i<<"\t"<<node_loc[0]<<"\t"<<node_loc[1]<<"\n";
        }
        node_file->close();
        std::string full_name = handler.GetTestOutputDirectory("")+"temp.";
        std::string command   = "./bin/triangle -e " + full_name + "node";
        system(command.c_str());
        
        TrianglesMeshReader<2,2> mesh_reader2(full_name+"1");
        ConformingTetrahedralMesh<2,2> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh2.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh2.GetNumBoundaryElements());
        
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh2.GetNumElements());
        
        //Test to see whether triangle/ tetgen is renumbering the nodes
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            const c_vector<double, 2> node_loc1 = mesh.GetNode(i)->rGetLocation();
            const c_vector<double, 2> node_loc2 = mesh2.GetNode(i)->rGetLocation();
            
            for (int j=0; j<2; j++)
            {
                TS_ASSERT_DELTA(node_loc1[j],node_loc2[j],1e-6);
            }
        }
    }
    
    void TestOperationOfTetgen() throw (Exception)
    {
        OutputFileHandler handler("");
        
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        double volume=mesh.CalculateMeshVolume();
        TS_ASSERT_DELTA(1, volume, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),375U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(),1626U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),390U);
        
        out_stream node_file=handler.OpenOutputFile("temp.node");
        (*node_file)<<mesh.GetNumNodes()<<"\t3\t0\t0\n";
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            const c_vector<double, 3> node_loc = mesh.GetNode(i)->rGetLocation();
            (*node_file)<<i<<"\t"<<node_loc[0]<<"\t"<<node_loc[1]<<"\t"<<node_loc[2]<<"\n";
        }
        node_file->close();
        std::string full_name = handler.GetTestOutputDirectory("")+"temp.";
        std::string command   = "./bin/tetgen -e " + full_name + "node";
        system(command.c_str());
        
        TrianglesMeshReader<3,3> mesh_reader2(full_name+"1");
        ConformingTetrahedralMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh2.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh2.GetNumBoundaryElements());
        
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh2.GetNumElements() + 2 );
        
        //Test to see whether triangle/ tetgen is renumbering the nodes
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            const c_vector<double, 3> node_loc1 = mesh.GetNode(i)->rGetLocation();
            const c_vector<double, 3> node_loc2 = mesh2.GetNode(i)->rGetLocation();
            
            for (int j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(node_loc1[j],node_loc2[j],1e-6);
            }
        }
    }
    
    // test 3d remesh - very similar test to TestOperationOfTetgenMoveNodes above, but
    // uses mesh.Remesh() instead of calling tetgen from here
    void TestRemesh3dMoveNodes() throw (Exception)
    {
        OutputFileHandler handler("");
        
        TrianglesMeshReader<3,3> mesh_reader2("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> old_mesh;
        old_mesh.ConstructFromMeshReader(mesh_reader2);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            Point<3> point = mesh.GetNode(i)->GetPoint();
            Point<3> old_mesh_point = old_mesh.GetNode(i)->GetPoint();
            for (int j=0; j<3; j++)
            {
                if (fabs(point[j]-0.0) >1e-6 && fabs(point[j]-1.0) >1e-6)
                {
                    point.SetCoordinate(j, point[j]+9e-2);
                    old_mesh_point.SetCoordinate(j, old_mesh_point[j]+9e-2);
                }
            }
            
            mesh.GetNode(i)->SetPoint(point);
            old_mesh.GetNode(i)->SetPoint(old_mesh_point);
        }
        mesh.RefreshMesh();
        old_mesh.RefreshMesh();

        double old_volume=mesh.CalculateMeshVolume();
        TS_ASSERT_DELTA(1, old_volume, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),375U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(),1626U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),390U);
        
        NodeMap map(mesh.GetNumNodes());
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), old_mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), old_mesh.GetNumBoundaryElements());
        
        TS_ASSERT_EQUALS(mesh.GetNumElements()+1, old_mesh.GetNumElements());
        
        //Test to see whether triangle/ tetgen is renumbering the nodes
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            // the map turns out to be the identity map in this test
            TS_ASSERT_EQUALS(map.GetNewIndex(i),i);
            
            const c_vector<double, 3> node_loc1 = mesh.GetNode(map.GetNewIndex(i))->rGetLocation();
            const c_vector<double, 3> node_loc2 = old_mesh.GetNode(i)->rGetLocation();
            
            for (int j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(node_loc1[j],node_loc2[j],1e-6);
            }
        }    
        
        double new_volume=mesh.CalculateMeshVolume();
        TS_ASSERT_DELTA(old_volume, new_volume, 1e-7);    
    }
    

    void TestOperationOfTetgenMoveNodes() throw (Exception)
    {
        OutputFileHandler handler("");
        
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            Point<3> point=mesh.GetNode(i)->GetPoint();
            for (int j=0; j<3; j++)
            {
                if (fabs(point[j]-0.0) >1e-6 && fabs(point[j]-1.0) >1e-6)
                {
                    point.SetCoordinate(j, point[j]+9e-2);
                }
            }
            
            mesh.GetNode(i)->SetPoint(point);
        }
        mesh.RefreshMesh();

        
        double volume=mesh.CalculateMeshVolume();
        TS_ASSERT_DELTA(1, volume, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),375U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(),1626U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),390U);
        
        out_stream node_file=handler.OpenOutputFile("temp.node");
        (*node_file)<<mesh.GetNumNodes()<<"\t3\t0\t0\n";
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            const c_vector<double, 3> node_loc=mesh.GetNode(i)->rGetLocation();
            (*node_file)<<i<<"\t"<<node_loc[0]<<"\t"<<node_loc[1]<<"\t"<<node_loc[2]<<"\n";
        }
        node_file->close();
        std::string full_name = handler.GetTestOutputDirectory("")+"temp.";
        std::string command   = "./bin/tetgen -e " + full_name + "node";
        system(command.c_str());
        
        TrianglesMeshReader<3,3> mesh_reader2(full_name+"1");
        ConformingTetrahedralMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh2.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh2.GetNumBoundaryElements());
        
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh2.GetNumElements() + 1);
        
        //Test to see whether triangle/ tetgen is renumbering the nodes
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            const c_vector<double, 3> node_loc1=mesh.GetNode(i)->rGetLocation();
            const c_vector<double, 3> node_loc2=mesh2.GetNode(i)->rGetLocation();
            
            for (int j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(node_loc1[j],node_loc2[j],1e-6);
            }
        }
    }


    void TestRemeshWithDeletions() throw (Exception)
    {
        OutputFileHandler handler("");
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        double area=mesh.CalculateMeshVolume();
        const int node_index=432;
        const int target_index=206;
        
        
        mesh.MoveMergeNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(area, mesh.CalculateMeshVolume(), 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),mesh.GetNumNodes()+1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());
        
        out_stream node_file=handler.OpenOutputFile("temp.node");
        (*node_file)<<mesh.GetNumNodes()<<"\t2\t0\t0\n";
        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
        {
            if (!mesh.GetNode(i)->IsDeleted())
            {
                const c_vector<double, 2> node_loc=mesh.GetNode(i)->rGetLocation();
                (*node_file)<<i<<"\t"<<node_loc[0]<<"\t"<<node_loc[1]<<"\n";
            }
        }
        
        unsigned new_index = 0;
        NodeMap map(mesh.GetNumAllNodes());
        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
        {
            if (mesh.GetNode(i)->IsDeleted())
            {
                map.SetDeleted(i);
            }
            else
            {
                map.SetNewIndex(i,new_index);
                new_index++;
            }
        }
        
        TS_ASSERT_EQUALS(new_index,mesh.GetNumNodes());
        
        node_file->close();
        std::string full_name = handler.GetTestOutputDirectory("")+"temp.";
        std::string command   = "./bin/triangle -e " + full_name + "node";
        system(command.c_str());
        
        TrianglesMeshReader<2,2> mesh_reader2(full_name+"1");
        ConformingTetrahedralMesh<2,2> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh2.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh2.GetNumBoundaryElements());
        
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh2.GetNumElements());
        
        //Test to see whether triangle/ tetgen is renumbering the nodes
        
        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
        {
            if (mesh.GetNode(i)->IsDeleted())
            {
                TS_ASSERT_THROWS_ANYTHING(map.GetNewIndex(i));
            }
            else
            {
                const c_vector<double, 2> node_loc1=mesh.GetNode(i)->rGetLocation();
                int another_new_index = map.GetNewIndex(i);
                const c_vector<double, 2> node_loc2=mesh2.GetNode(another_new_index)->rGetLocation();
                
                for (int j=0; j<2; j++)
                {
                    TS_ASSERT_DELTA(node_loc1[j],node_loc2[j],1e-6);
                }
            }
        }
    }
    
    void TestRemeshWithMethod2D() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        double area=mesh.CalculateMeshVolume();
        const int node_index=432;
        const int target_index=206;
        
        unsigned num_nodes_before=mesh.GetNumNodes();
        unsigned num_elements_before=mesh.GetNumElements();
        unsigned num_boundary_elements_before=mesh.GetNumBoundaryElements();
        
        mesh.MoveMergeNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(area, mesh.CalculateMeshVolume(), 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),mesh.GetNumNodes()+1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());
        
        NodeMap map(1);
        mesh.ReMesh(map);
        
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());
        
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), num_elements_before-2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), num_nodes_before-1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), num_boundary_elements_before);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(),area,1e-6);
    }
    
    void TestRemeshCrinklyNonVoronoi() throw (Exception)
    {
        OutputFileHandler handler("");
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/crinkly");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
//        double area=mesh.CalculateMeshVolume();
//
//        int num_nodes_before=mesh.GetNumNodes();
//        int num_elements_before=mesh.GetNumElements();
//        int num_boundary_elements_before=mesh.GetNumBoundaryElements();
//
        NodeMap map(1);
        mesh.ReMesh(map);
        
//   	    TS_ASSERT_EQUALS(mesh.GetNumAllElements(), num_elements_before);
//        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),num_nodes_before);
//        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), num_boundary_elements_before);
//   		TS_ASSERT_DELTA(mesh.CalculateMeshVolume(),area,1e-6);
//
    }
    
    void TestCreateMirrorCells() throw (Exception)
    {   // note that elements are not created (and boundary elements are not changed)
        // this just creates a set of new nodes.
    	unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        
        double x0 = 0.0;
        double x1 = crypt_width;
        std::vector<std::vector<unsigned> > image_map;
        
        // Test that the method throws if the mesh has nodes outside {x0,x1}
        TS_ASSERT_THROWS_ANYTHING(image_map = p_mesh->CreateMirrorNodes(x0,x1));
        
        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator2(cells_across, cells_up, 0u, true);
        p_mesh=generator2.GetMesh();
        
        TS_ASSERT_THROWS_NOTHING(image_map = p_mesh->CreateMirrorNodes(x0,x1));
        
        //Output2DMeshToFile(p_mesh, "node_positions.dat");
        
        std::vector<unsigned> left_original = image_map[0];
        std::vector<unsigned> left_images = image_map[1];
        std::vector<unsigned> right_original = image_map[2];
        std::vector<unsigned> right_images = image_map[3];
        
        //std::cout << "Left size = " << left_original.size()<< ", Right size = " << right_original.size()<< "\n" << std::flush;
        // Check the vectors are the right size...
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
    }
    
    void TestReconstructCylindricalMesh() throw (Exception)
    {   // this takes in a new mesh created using the mirror function above
        // and a ReMesh call, then removes nodes, elements and boundary elements.
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 0;
        
        double xLeft = 0.0;
        double xRight = crypt_width;
        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        
        // Create a mirrored load of nodes for the normal remesher to work with.
        std::vector<std::vector<unsigned> > image_map = p_mesh->CreateMirrorNodes(xLeft, xRight);
    
        std::vector<unsigned> left_original = image_map[0];
        std::vector<unsigned> left_images = image_map[1];
        std::vector<unsigned> right_original = image_map[2];
        std::vector<unsigned> right_images = image_map[3];
        
        // Call the normal re-mesh
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);
        
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
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        
        CancerParameters *p_params = CancerParameters::Instance();
        double x0 = 0.0;
        double x1 = p_params->GetCryptWidth();
        p_mesh->CylindricalReMesh(x0,x1);

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
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        
        CancerParameters *p_params = CancerParameters::Instance();
        
        double x0 = 0.0;
        double x1 = p_params->GetCryptWidth();
        p_mesh->CylindricalReMesh(x0,x1);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(),2*cells_across);

        //Output2DMeshToFile(p_mesh, "node_positions.dat");
    }
    
    void TestGetDistanceBetweenCyclindricalPoints() throw (Exception)
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;
        
        // Set up a mesh which can be mirrored (no ghosts in this case)
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        
        CancerParameters *p_params = CancerParameters::Instance();
        
        double x0 = 0.0;
        double x1 = p_params->GetCryptWidth();

        c_vector<double, 2> location1 = p_mesh->GetNode(1)->rGetLocation();
        c_vector<double, 2> location2 = p_mesh->GetNode(4)->rGetLocation();

        // test a normal distance calculation
        double distance = p_mesh->GetDistanceBetweenCylindricalPoints(location1, location2, x0,x1);
        TS_ASSERT_DELTA(distance, 1.0, 1e-7);
        
        // test a periodic calculation
        location1[0] = 0.5;
        location1[1] = 3.0;
        location2[0] = 2.5;
        location2[1] = 4.0;
        distance = p_mesh->GetDistanceBetweenCylindricalPoints(location1, location2, x0,x1);
        TS_ASSERT_DELTA(distance, sqrt(2.0), 1e-7);
        
        // test a periodic calculation where points need to be swapped
        location1[0] = 2.5;
        location1[1] = 4.0;
        location2[0] = 0.5;
        location2[1] = 3.0;
        distance = p_mesh->GetDistanceBetweenCylindricalPoints(location1, location2, x0,x1);
        TS_ASSERT_DELTA(distance, sqrt(2.0), 1e-7);
    }
    

    
};

#endif /*TESTREMESH_HPP_*/
