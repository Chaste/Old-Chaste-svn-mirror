#ifndef TESTPOINTINCLUSION_HPP_
#define TESTPOINTINCLUSION_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>

class TestPointInclusion : public CxxTest::TestSuite
{
public:

	void TestPointInElement1D()
	{
        std::vector<Node<1>*> nodes1d;
        nodes1d.push_back(new Node<1>(0, false, 2.0));
        nodes1d.push_back(new Node<1>(1, false, 2.5));
        Element<1,1> element1d(INDEX_IS_NOT_USED, nodes1d);
        
        Point<1> in_point(2.25);
        Point<1> on_point(2.00);
        Point<1> out_point(1.25);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(in_point), true);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(out_point), false);
        
        delete nodes1d[0];
        delete nodes1d[1];
	}
    
    
   void TestPointinMesh1D(void)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<1> point1(0.15);
        Point<1> point2(-0.1);
        Point<1> point3(0.2);
        TS_ASSERT_EQUALS(mesh.GetContainingElement(point1),1U);
        TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElement(point2));
        TS_ASSERT_EQUALS(mesh.GetContainingElement(point3),1U);  //in elements 1 and 2
    }
  
        
   	void TestPointInElement2D()
	{
        std::vector<Node<2>*> nodes2d;
        nodes2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes2d.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes2d.push_back(new Node<2>(2, false, 0.0, 3.0));
        Element<2,2> element2d(INDEX_IS_NOT_USED, nodes2d);
        //c_matrix<double,2,2> det=*(element2d.GetJacobian());
        //std::cout<<det<<"\n";
        delete nodes2d[0];
        delete nodes2d[1];
        delete nodes2d[2];
        
         	
	}
    
    

  
        
   	void TestPointInElement3D()
	{
        std::vector<Node<3>*> nodes3d;
        nodes3d.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes3d.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        delete nodes3d[0];
        delete nodes3d[1];
        delete nodes3d[2];
        delete nodes3d[3];
	}        
};

#endif /*TESTPOINTINCLUSION_HPP_*/
