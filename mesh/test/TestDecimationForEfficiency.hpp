#ifndef _TESTDECIMATIONFOREFFICIENCY_HPP_
#define _TESTDECIMATIONFOREFFICIENCY_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include "Decimator.hpp"
#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>
#include <fstream>

class TestDecimationForEfficiency : public CxxTest::TestSuite
{
public:

    void TestDecimateHeart() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/heart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59799, 1.0e-5);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 63885U);
        
        Decimator<3> decimator;
        
        decimator.Initialise(&mesh);
        TS_ASSERT_DELTA(decimator.GetVolumeLeakage(), 1e-5, 1e-20)
        //decimator.SetVolumeLeakage(1e-4);
        //decimator.SetVolumeLeakage(1e-3);
        //decimator.SetVolumeLeakage(1e-2);
        decimator.SetVolumeLeakage(5e-2);
        //decimator.SetVolumeLeakage(1e-1);
        
        
        decimator.Decimate();
        
        
 	    //std::cout<<mesh.CalculateMeshVolume()<<"\n";
        //TS_ASSERT_EQUALS(mesh.GetNumNodes(), 19516U);//VolumeLeakage(1e-4);  
        //TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59786, 1.0e-5);//VolumeLeakage(1e-4);
        //TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10555U);//VolumeLeakage(1e-3);
        //TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.59321, 1.0e-4);//VolumeLeakage(1e-3);
        //TS_ASSERT_EQUALS(mesh.GetNumNodes(), 1545U);//VolumeLeakage(1e-2);
        //TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.53256, 1.0e-4);//VolumeLeakage(1e-2);
 	    TS_ASSERT_EQUALS(mesh.GetNumNodes(), 173U);//VolumeLeakage(5e-2);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.40447, 1.0e-4);//VolumeLeakage(5e-2);
 	    //TS_ASSERT_EQUALS(mesh.GetNumNodes(), 62U);//VolumeLeakage(1e-1);
        //TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 6.67374, 1.0e-4);//VolumeLeakage(1e-1);
 	
        TrianglesMeshWriter<3,3> mesh_writer("", "HeartDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
    
      
        
        
        
   
    
   
 

    
};
#endif //_TESTDECIMATIONFOREFFICIENCY_HPP_
