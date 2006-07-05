#ifndef _TESTNONLINEARELASTICPDE_HPP_
#define _TESTNONLINEARELASTICPDE_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"

#include "CompressibleIsotropicMooneyRivlinMaterial.hpp"
#include "NonlinearElasticityPde.hpp"

class TestNonlinearElasticityPde : public CxxTest::TestSuite 
{
public:
	void testNonlinearElasticityPde2D()
	{
		TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		double c1 = 12.23;
		CompressibleIsotropicMooneyRivlinMaterial<2> material(c1);
		material.SetDensity(1.5);
	
		c_vector<double, 2> gravity;
		gravity(0) = 1;
		gravity(1) = 2;
		NonlinearElasticityPde<2,2> pde(gravity);

		for (int i=0; i<mesh.GetNumElements(); i++)		 
		{
			mesh.SetMaterialToElement(i, &material);		
		}

//		for (int i=0; i<mesh.GetNumElements(); i++)		 
//		{
//			MatrixDouble I = MatrixDouble::Identity(2);						 
//			MatrixDouble T = pde.ComputeStress( mesh.GetElement(i), I);
//
//			c_vector<double, SPACE_DIM> rhoG = pde.ComputeGravityForceTerm(  mesh.GetElement(i) );
//			
//			TS_ASSERT_EQUALS( T.Rows()   ,2);
//			TS_ASSERT_EQUALS( T.Columns(),2);
//			TS_ASSERT_EQUALS( rhoG.Size(),2);
//			
//			TS_ASSERT_DELTA( T(0,0),  0, 1e-12);
//			TS_ASSERT_DELTA( T(1,0),  0, 1e-12);
//			TS_ASSERT_DELTA( T(0,1),  0, 1e-12);
//			TS_ASSERT_DELTA( T(1,1),  0, 1e-12);
//			TS_ASSERT_DELTA( rhoG(0), 1.5, 1e-12);							
//			TS_ASSERT_DELTA( rhoG(1), 3.0, 1e-12);							
//		}
	}

	void testNonlinearElasticityPde3D()
	{
		TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
		ConformingTetrahedralMesh<3,3> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		double c1 = 12.23;
		double c2 = 13.34;
		CompressibleIsotropicMooneyRivlinMaterial<3> material(c1,c2);
		material.SetDensity(1.5);
	
		c_vector<double, 3> gravity;
		gravity(0) = 1;
		gravity(1) = 2;
		gravity(2) = 3;
		NonlinearElasticityPde<3,3> pde(gravity);

		for (int i=0; i<mesh.GetNumElements(); i++)		 
		{
			mesh.SetMaterialToElement(i, &material);		
		}

//		for (int i=0; i<mesh.GetNumElements(); i++)		 
//		{
//			MatrixDouble I = MatrixDouble::Identity(3);						 
//			MatrixDouble T = pde.ComputeStress( mesh.GetElement(i), I);
//
//			c_vector<double, SPACE_DIM> rhoG = pde.ComputeGravityForceTerm(  mesh.GetElement(i) );
//			
//			TS_ASSERT_EQUALS( T.Rows()   ,3);
//			TS_ASSERT_EQUALS( T.Columns(),3);
//			TS_ASSERT_EQUALS( rhoG.Size(),3);
//			
//			TS_ASSERT_DELTA( T(0,0),  0, 1e-12);
//			TS_ASSERT_DELTA( T(0,1),  0, 1e-12);
//			TS_ASSERT_DELTA( T(0,2),  0, 1e-12);
//			TS_ASSERT_DELTA( T(1,0),  0, 1e-12);
//			TS_ASSERT_DELTA( T(1,1),  0, 1e-12);
//			TS_ASSERT_DELTA( T(1,2),  0, 1e-12);
//			TS_ASSERT_DELTA( T(2,0),  0, 1e-12);
//			TS_ASSERT_DELTA( T(2,1),  0, 1e-12);
//			TS_ASSERT_DELTA( T(2,2),  0, 1e-12);
//			TS_ASSERT_DELTA( rhoG(0), 1.5, 1e-12);							
//			TS_ASSERT_DELTA( rhoG(1), 3.0, 1e-12);							
//			TS_ASSERT_DELTA( rhoG(2), 4.5, 1e-12);							
//		}
	}
};

#endif //_TESTNONLINEARELASTICPDE_HPP_
