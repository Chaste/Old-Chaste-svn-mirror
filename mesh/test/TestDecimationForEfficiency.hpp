#ifndef _TESTDECIMATIONFOREFFICIENCY_HPP_
#define _TESTDECIMATIONFOREFFICIENCY_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include "RandomDecimator.hpp"
#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>
#include <fstream>

class TestDecimationForEfficiency : public CxxTest::TestSuite
{
public:

    void TestRandomDecimator3D() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/heart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59799, 1.0e-5);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 63885U);
        
        RandomNumberGenerator::Instance();
        RandomDecimator<3> decimator;
        decimator.Initialise(&mesh);
        decimator.SetVolumeLeakage(1e-5);
        
        
        decimator.Decimate();
        RandomNumberGenerator::Destroy();
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21894U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59799, 1.0e-5);
        TrianglesMeshWriter<3,3> mesh_writer("", "HeartDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
        
        
    }
};
#endif //_TESTDECIMATIONFOREFFICIENCY_HPP_
