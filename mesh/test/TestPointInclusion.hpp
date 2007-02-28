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

    void TestPointWeightsInElement1D()
    {
        std::vector<Node<1>*> nodes1d;
        nodes1d.push_back(new Node<1>(0, false, 2.0));
        nodes1d.push_back(new Node<1>(1, false, 2.5));
        Element<1,1> element1d(INDEX_IS_NOT_USED, nodes1d);
        
        Point<1> in_point(2.25);
        Point<1> on_point(2.00);
        Point<1> out_point(1.25);
        c_vector <double, 2> weights;
        weights=element1d.CalculateInterpolationWeights(on_point);
        TS_ASSERT_EQUALS(weights[0], 1.0);
        TS_ASSERT_EQUALS(weights[1], 0.0);
        weights=element1d.CalculateInterpolationWeights(in_point);
        TS_ASSERT_EQUALS(weights[0], 0.5);
        TS_ASSERT_EQUALS(weights[1], 0.5);
        weights=element1d.CalculateInterpolationWeights(out_point);
        //1.25 = 2.5*2 -1.5 * 2.5
        TS_ASSERT_EQUALS(weights[0], 2.5);
        TS_ASSERT_EQUALS(weights[1], -1.5);
        
      
        
        delete nodes1d[0];
        delete nodes1d[1];
    }
    
    void TestPointInElement1D()
    {
        std::vector<Node<1>*> nodes1d;
        nodes1d.push_back(new Node<1>(0, false, 2.0));
        nodes1d.push_back(new Node<1>(1, false, 2.5));
        Element<1,1> element1d(INDEX_IS_NOT_USED, nodes1d);
        
        Point<1> in_point(2.25);
        Point<1> on_point(2.00);
        Point<1> out_point(1.25);
        bool strict=true;
        TS_ASSERT_EQUALS(element1d.IncludesPoint(in_point), true);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(on_point, strict), false);
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
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1),1U);
        TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(point2));
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point3),1U);  //in elements 1 and 2
        
        std::vector<unsigned> indices;
        indices=mesh.GetContainingElementIndices(point1);
        TS_ASSERT_EQUALS(indices.size(), 1U);
        TS_ASSERT_EQUALS(indices[0], 1U);

        indices=mesh.GetContainingElementIndices(point2);
        TS_ASSERT_EQUALS(indices.size(), 0U);

        indices=mesh.GetContainingElementIndices(point3);
        TS_ASSERT_EQUALS(indices.size(), 2U);
        TS_ASSERT_EQUALS(indices[0], 1U);
        TS_ASSERT_EQUALS(indices[1], 2U);
    }
  
        
   	void TestPointWeightsAndInclusion2D()
	{
        std::vector<Node<2>*> nodes2d;
        nodes2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes2d.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes2d.push_back(new Node<2>(2, false, 0.0, 3.0));
        Element<2,2> element2d(INDEX_IS_NOT_USED, nodes2d);
       
        Point<2> on_point(0., 2.);
        c_vector <double, 3> weights;
        bool strict=true;
        TS_ASSERT_EQUALS(element2d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element2d.IncludesPoint(on_point, strict), false);
        weights=element2d.CalculateInterpolationWeights(on_point);
        TS_ASSERT_DELTA(weights[0], 1.0/3.0, 1e-5);
        TS_ASSERT_DELTA(weights[1], 0.0, 1e-5);
        TS_ASSERT_DELTA(weights[2], 2.0/3.0, 1e-5);

        Point<2> in_point(1., 1.);
        TS_ASSERT_EQUALS(element2d.IncludesPoint(in_point), true);
        weights=element2d.CalculateInterpolationWeights(in_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(0.0, weights[2]);

        Point<2> out_point(1., 0.);
        TS_ASSERT_EQUALS(element2d.IncludesPoint(out_point), false);
        weights=element2d.CalculateInterpolationWeights(out_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(weights[2], 0.0);

        delete nodes2d[0];
        delete nodes2d[1];
        delete nodes2d[2];
        
         	
	}
    
    
    void TestPointinMesh2D(void)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<2> point1(0.051, 0.051);
        Point<2> point2(0.2,0.2);
        Point<2> point3(0.05, 0.05); //Node 60 of mesh
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1),110U);
        TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(point2));
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point3),89U);  //in elements 89,90,91,108,109, 110
        
        
        std::vector<unsigned> indices;
        indices=mesh.GetContainingElementIndices(point1);
        TS_ASSERT_EQUALS(indices.size(), 1U);
        TS_ASSERT_EQUALS(indices[0], 110U);

        indices=mesh.GetContainingElementIndices(point2);
        TS_ASSERT_EQUALS(indices.size(), 0U);

        indices=mesh.GetContainingElementIndices(point3);
        TS_ASSERT_EQUALS(indices.size(), 6U);
        TS_ASSERT_EQUALS(indices[0], 89U);
        TS_ASSERT_EQUALS(indices[1], 90U);
        TS_ASSERT_EQUALS(indices[5], 110U);
    }
  
  
        
   	void TestPointInElement3D()
	{
        std::vector<Node<3>*> nodes3d;
        nodes3d.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes3d.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element3d(INDEX_IS_NOT_USED, nodes3d);
    
        bool strict=true;
        Point<3> on_point(0., 0.2, 0.);
        c_vector <double, 4> weights;
        TS_ASSERT_EQUALS(element3d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(on_point, strict), false);
        weights=element3d.CalculateInterpolationWeights(on_point);
        TS_ASSERT_DELTA(weights[0], 0.8, 1e-5);
        TS_ASSERT_DELTA(weights[1], 0.0, 1e-5);
        TS_ASSERT_DELTA(weights[2], 0.2, 1e-5);
        TS_ASSERT_DELTA(weights[3], 0.0, 1e-5);

        Point<3> in_point(.25, .25, .25);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(in_point), true);
        weights=element3d.CalculateInterpolationWeights(in_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(0.0, weights[2]);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);

        Point<3> out_point(0.1, -10., 0.1);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(out_point), false);
        weights=element3d.CalculateInterpolationWeights(out_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(weights[2], 0.0);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);
 
 
 
        delete nodes3d[0];
        delete nodes3d[1];
        delete nodes3d[2];
        delete nodes3d[3];
	}
    
    void TestPointinMesh3D(void)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<3> point1(0.051, 0.051,0.051);
        Point<3> point2(0.2,0.2,0.2);
        Point<3> point3(0.050000000000000003,  0.050000000000000003,  0.050000000000000003);
        //Node 665 of mesh
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1),2992U);
        TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(point2));
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point3),2044U);  
        /*in elements 2044, 2047. 2058, 2192, 2268, 2286, 2392, 2414, 2415,
         * 2424, 2426, 2452, 2661, 2704, 2734, 2745, 2846, 2968, 2990, 2992,
         * 3015, 3022, 3024, 3026
         */
         
        //should throw because vertex is not strictly contained in any element
        TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(point3, true));  
        
        
        std::vector<unsigned> indices;
        indices=mesh.GetContainingElementIndices(point1);
        TS_ASSERT_EQUALS(indices.size(), 1U);
        TS_ASSERT_EQUALS(indices[0], 2992U);

        indices=mesh.GetContainingElementIndices(point2);
        TS_ASSERT_EQUALS(indices.size(), 0U);

        indices=mesh.GetContainingElementIndices(point3);
        TS_ASSERT_EQUALS(indices.size(), 24U);
        TS_ASSERT_EQUALS(indices[0], 2044U);
        TS_ASSERT_EQUALS(indices[1], 2047U);
        TS_ASSERT_EQUALS(indices[5], 2286U);   
        TS_ASSERT_EQUALS(indices[23], 3026U); 
         
    }
          
};

#endif /*TESTPOINTINCLUSION_HPP_*/
