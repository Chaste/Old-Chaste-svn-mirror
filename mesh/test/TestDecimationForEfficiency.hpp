#ifndef _TESTDECIMATIONFOREFFICIENCY_HPP_
#define _TESTDECIMATIONFOREFFICIENCY_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include "RandomDecimator.hpp"
#include "SequenceDecimator.hpp"
#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>
#include <fstream>


/*
template <unsigned SPACE_DIM>
class FixedNodeDecimator : public Decimator<SPACE_DIM>
{
private:
	unsigned mEndNumberNodes;
	
protected:	
   	bool ExtraStoppingCondition()
	{

		return (this->mQueue.size() == mEndNumberNodes);
	}

public:
	FixedNodeDecimator(unsigned endNumberNodes)
	{
		mEndNumberNodes=endNumberNodes;
	}

};
*/

class TestDecimationForEfficiency : public CxxTest::TestSuite
{
public:

    void xTestDecimateHeart() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/heart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59799, 1.0e-5);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 63885U);
        
        RandomNumberGenerator::Instance();
        Decimator<3> decimator;
        
        decimator.Initialise(&mesh);
        TS_ASSERT_DELTA(decimator.GetVolumeLeakage(), 1e-5, 1e-20)
        decimator.SetVolumeLeakage(1e-4);
        
        
        //decimator.Decimate(true);
        decimator.Decimate();
        RandomNumberGenerator::Destroy();
        
        //TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21894U); //Random
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 19516U);
        //TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21800U); //Sequence
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59786, 1.0e-5);
        TrianglesMeshWriter<3,3> mesh_writer("", "HeartDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
    
    void xTestAgressiveDecimateHeart() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/heart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59799, 1.0e-5);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 63885U);
        
        RandomNumberGenerator::Instance();
        Decimator<3> decimator;
        
        decimator.Initialise(&mesh);
        TS_ASSERT_DELTA(decimator.GetVolumeLeakage(), 1e-5, 1e-20)
        decimator.SetVolumeLeakage(1e-3);
        
        
        decimator.Decimate();
        //decimator.Decimate(true);
        RandomNumberGenerator::Destroy();
        
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10555U);
        
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59321, 1.0e-4);
        TrianglesMeshWriter<3,3> mesh_writer("", "HeartAgressiveDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
    
    
    void TestMostAgressiveDecimateHeart() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/heart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59799, 1.0e-5);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 63885U);
        
        RandomNumberGenerator::Instance();
        Decimator<3> decimator;
        
        decimator.Initialise(&mesh);
        TS_ASSERT_DELTA(decimator.GetVolumeLeakage(), 1e-5, 1e-20)
        decimator.SetVolumeLeakage(1e-2);
        
        
        //decimator.Decimate(true);
        decimator.Decimate();
        RandomNumberGenerator::Destroy();
        
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 1545U);
        
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.53256, 1.0e-4);
        TrianglesMeshWriter<3,3> mesh_writer("", "HeartMostAgressiveDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
  
 
 

    
};
#endif //_TESTDECIMATIONFOREFFICIENCY_HPP_
