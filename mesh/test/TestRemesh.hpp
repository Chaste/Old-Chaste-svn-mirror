#ifndef TESTREMESH_HPP_
#define TESTREMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "OutputFileHandler.hpp"
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
        		
	}
};

#endif /*TESTREMESH_HPP_*/
