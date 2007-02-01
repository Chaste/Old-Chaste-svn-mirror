#ifndef _TESTDECIMATIONFOREFFICIENCY_HPP_
#define _TESTDECIMATIONFOREFFICIENCY_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "RandomDecimator.hpp"
#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>
#include <fstream>

class TestDecimationForEfficiency : public CxxTest::TestSuite
{
public:
    void xTestRandomDecimator2D() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        RandomDecimator<2> decimator;
        decimator.Initialise(&mesh);
        
        decimator.DecimateAnimate("RandomAnimation");
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
    }
    void TestRandomDecimator3D() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/heart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59799, 1.0e-5);
        Decimator<3> decimator;
        decimator.Initialise(&mesh);
        
        decimator.Decimate();
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59799, 1.0e-5);
    }
};
#endif //_TESTDECIMATIONFOREFFICIENCY_HPP_
