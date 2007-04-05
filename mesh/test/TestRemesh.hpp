#ifndef TESTREMESH_HPP_
#define TESTREMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "OutputFileHandler.hpp"
#include "NodeMap.hpp"
#include <iostream>
#include <cmath>

#include <vector>

class TestRemesh : public CxxTest::TestSuite
{
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


    void TestRemeshWithDeletions()
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
    
    void TestRemeshWithMethod2D()
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
    
    void TestRemeshCrinklyNonVoronoi()
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
    
    
};

#endif /*TESTREMESH_HPP_*/
