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
   		TS_ASSERT_EQUALS(mesh.GetNumNodes(),77);
   		TS_ASSERT_EQUALS(mesh.GetNumElements(),141);
   		TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),11);
   		
   		out_stream node_file=handler.OpenOutputFile("temp.node");
   		(*node_file)<<mesh.GetNumNodes()<<"\t2\t0\t0\n";
   		for (int i=0; i<mesh.GetNumNodes(); i++)
   		{
   			Point<2> point=mesh.GetNodeAt(i)->rGetPoint();
   			(*node_file)<<i<<"\t"<<point[0]<<"\t"<<point[1]<<"\n";
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
 		
 		for (int i=0; i<mesh.GetNumNodes(); i++)
   		{
        	Point<2> point1=mesh.GetNodeAt(i)->rGetPoint();	
        	Point<2> point2=mesh2.GetNodeAt(i)->rGetPoint();
        	
        	for (int j=0; j<2; j++)
        	{
        		TS_ASSERT_DELTA(point1[j],point2[j],1e-6);
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
   		TS_ASSERT_EQUALS(mesh.GetNumNodes(),375);
   		TS_ASSERT_EQUALS(mesh.GetNumElements(),1626);
   		TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),390);
   		
   		out_stream node_file=handler.OpenOutputFile("temp.node");
   		(*node_file)<<mesh.GetNumNodes()<<"\t3\t0\t0\n";
   		for (int i=0; i<mesh.GetNumNodes(); i++)
   		{
   			Point<3> point=mesh.GetNodeAt(i)->rGetPoint();
   			(*node_file)<<i<<"\t"<<point[0]<<"\t"<<point[1]<<"\t"<<point[2]<<"\n";
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
 	
 		TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh2.GetNumElements());
 		
 		//Test to see whether triangle/ tetgen is renumbering the nodes
 		
 		for (int i=0; i<mesh.GetNumNodes(); i++)
   		{
        	Point<3> point1=mesh.GetNodeAt(i)->rGetPoint();	
        	Point<3> point2=mesh2.GetNodeAt(i)->rGetPoint();
        	
        	for (int j=0; j<3; j++)
        	{
        		TS_ASSERT_DELTA(point1[j],point2[j],1e-6);
        	}
   		}
        		
	}
	
	void TestOperationOfTetgenMoveNodes() throw (Exception)
	{
		OutputFileHandler handler("");
		
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        for (int i=0; i<mesh.GetNumNodes(); i++)
   		{
   			Point<3> point=mesh.GetNodeAt(i)->GetPoint();
   			for (int j=0; j<3; j++)
   			{
   				if (fabs(point[j]-0.0) >1e-6 && fabs(point[j]-1.0) >1e-6)
   				{
   					point.SetCoordinate(j, point[j]+9e-2);
   				}
   			}
   			
   			mesh.GetNodeAt(i)->SetPoint(point);
   		}
        mesh.RefreshMesh();
        
        
        
        double volume=mesh.CalculateMeshVolume();
   		TS_ASSERT_DELTA(1, volume, 1e-7);
   		TS_ASSERT_EQUALS(mesh.GetNumNodes(),375);
   		TS_ASSERT_EQUALS(mesh.GetNumElements(),1626);
   		TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),390);
   		
   		out_stream node_file=handler.OpenOutputFile("temp.node");
   		(*node_file)<<mesh.GetNumNodes()<<"\t3\t0\t0\n";
   		for (int i=0; i<mesh.GetNumNodes(); i++)
   		{
   			Point<3> point=mesh.GetNodeAt(i)->rGetPoint();
   			(*node_file)<<i<<"\t"<<point[0]<<"\t"<<point[1]<<"\t"<<point[2]<<"\n";
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
 	
 		TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh2.GetNumElements());
 		
 		//Test to see whether triangle/ tetgen is renumbering the nodes
 		
 		for (int i=0; i<mesh.GetNumNodes(); i++)
   		{
        	Point<3> point1=mesh.GetNodeAt(i)->rGetPoint();	
        	Point<3> point2=mesh2.GetNodeAt(i)->rGetPoint();
        	
        	for (int j=0; j<3; j++)
        	{
        		TS_ASSERT_DELTA(point1[j],point2[j],1e-6);
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
       
        
        mesh.SetNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(area, mesh.CalculateMeshVolume(), 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),mesh.GetNumNodes()+1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());
        
        out_stream node_file=handler.OpenOutputFile("temp.node");
   		(*node_file)<<mesh.GetNumNodes()<<"\t2\t0\t0\n";
   		for (int i=0; i<mesh.GetNumAllNodes(); i++)
   		{
   			if (!mesh.GetNodeAt(i)->IsDeleted())
   			{
   				Point<2> point=mesh.GetNodeAt(i)->rGetPoint();
   				(*node_file)<<i<<"\t"<<point[0]<<"\t"<<point[1]<<"\n";
   			}   			
   		}
   		
   		int new_index = 0;
   		NodeMap map(mesh.GetNumAllNodes());
   		for (int i=0; i<mesh.GetNumAllNodes(); i++)
   		{
   			if (mesh.GetNodeAt(i)->IsDeleted())
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
 		
 		for (int i=0; i<mesh.GetNumAllNodes(); i++)
   		{
   			if(mesh.GetNodeAt(i)->IsDeleted())
   			{
   				TS_ASSERT_THROWS_ANYTHING(map.GetNewIndex(i));
   			}
   			else
   			{
        		Point<2> point1=mesh.GetNodeAt(i)->rGetPoint();
        		int another_new_index = map.GetNewIndex(i);
        		Point<2> point2=mesh2.GetNodeAt(another_new_index)->rGetPoint();
        	
      		  	for (int j=0; j<2; j++)
        		{
        			TS_ASSERT_DELTA(point1[j],point2[j],1e-6);
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
       
        int num_nodes_before=mesh.GetNumNodes();
        
        mesh.SetNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(area, mesh.CalculateMeshVolume(), 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),mesh.GetNumNodes()+1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());
        
        NodeMap map(1);
        mesh.ReMesh(map);

		TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());

		TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 0);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),num_nodes_before-1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), 0);
   		
	}
        
	
	
	
};

#endif /*TESTREMESH_HPP_*/
