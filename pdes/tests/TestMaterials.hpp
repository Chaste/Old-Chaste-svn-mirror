#ifndef _TESTMATERIALS_HPP_
#define _TESTMATERIALS_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractMaterial.hpp"
#include "CompressibleIsotropicMooneyRivlinMaterial.hpp"
 
class TestMaterials : public CxxTest::TestSuite 
{
public:
	void testCompressibleIsotropicMooneyRivlinMaterial()
	{
		CompressibleIsotropicMooneyRivlinMaterial<1> material1(2);
		material1.SetDensity(1.4);
		TS_ASSERT_DELTA( material1.GetDensity(), 1.4, 1e-12);
		
		MatrixDouble G(1,1);
		G(0,0) = 2;

		MatrixDouble T1 = material1.ComputeStress(G);
		
		TS_ASSERT_EQUALS( T1.Rows(), 1);
		TS_ASSERT_EQUALS( T1.Columns(), 1);
		TS_ASSERT_DELTA(  T1(0,0), 4, 1e-4); 
		
		CompressibleIsotropicMooneyRivlinMaterial<2> material2(2);
		material2.SetDensity(1.5);
		TS_ASSERT_DELTA( material2.GetDensity(), 1.5, 1e-12);
		
		MatrixDouble H(2,2);
		H(0,0) = 1;
		H(0,1) = 2;
		H(1,0) = 3;
		H(1,1) = 4;
		
		MatrixDouble T2 = material2.ComputeStress(H);
		
		TS_ASSERT_EQUALS( T2.Rows(), 2);
		TS_ASSERT_EQUALS( T2.Columns(), 2);

		TS_ASSERT_DELTA( T2(0,0), -76, 1e-4); 
		TS_ASSERT_DELTA( T2(0,1),  56, 1e-4); 
		TS_ASSERT_DELTA( T2(1,0),  56, 1e-4); 
		TS_ASSERT_DELTA( T2(1,1), -36, 1e-4);  
		
				
		CompressibleIsotropicMooneyRivlinMaterial<3> material3(2,1);
		
		material3.SetDensity(1.3);
		TS_ASSERT_DELTA( material3.GetDensity(), 1.3, 1e-12);
		
		MatrixDouble F(3,3);
		F(0,0) = 1;
		F(0,1) = 2;
		F(0,2) = 3;  
		F(1,0) = 4; 
		F(1,1) = 5;
		F(1,2) = 6;
		F(2,0) = 7; 
		F(2,1) = 8;
		F(2,2) = 8; //not 9
		
		MatrixDouble T = material3.ComputeStress(F);
		
		TS_ASSERT_EQUALS( T.Rows(), 3);
		TS_ASSERT_EQUALS( T.Columns(), 3);

		TS_ASSERT_DELTA( T(0,0), -688, 1e-4); 
		TS_ASSERT_DELTA( T(0,1), 1460, 1e-4); 
		TS_ASSERT_DELTA( T(0,2), -814, 1e-4); 
		TS_ASSERT_DELTA( T(1,0), 1460, 1e-4);  
		TS_ASSERT_DELTA( T(1,1),-2086, 1e-4); 
		TS_ASSERT_DELTA( T(1,2),  808, 1e-4); 
		TS_ASSERT_DELTA( T(2,0), -814, 1e-4); 
		TS_ASSERT_DELTA( T(2,1),  808, 1e-4); 
		TS_ASSERT_DELTA( T(2,2), -110, 1e-4); 
		 		   
 		MatrixDouble I = MatrixDouble::Identity(3);
 		T = material3.ComputeStress(I);

		TS_ASSERT_DELTA( T(0,0),  0, 1e-4); 
		TS_ASSERT_DELTA( T(0,1),  0, 1e-4); 
		TS_ASSERT_DELTA( T(0,2),  0, 1e-4); 
		TS_ASSERT_DELTA( T(1,0),  0, 1e-4);  
		TS_ASSERT_DELTA( T(1,1),  0, 1e-4); 
		TS_ASSERT_DELTA( T(1,2),  0, 1e-4); 
		TS_ASSERT_DELTA( T(2,0),  0, 1e-4); 
		TS_ASSERT_DELTA( T(2,1),  0, 1e-4); 
		TS_ASSERT_DELTA( T(2,2),  0, 1e-4); 
		
	}
};


#endif //_TESTMATERIALS_HPP_
