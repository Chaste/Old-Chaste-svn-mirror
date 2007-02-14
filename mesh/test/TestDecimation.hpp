#ifndef _TESTDECIMATION_HPP_
#define _TESTDECIMATION_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include "Decimator.hpp"
#include "QualityDecimator.hpp"
#include "SequenceDecimator.hpp"
#include "MinimumElementDecimator.hpp"
#include "LinearFunctionDecimator.hpp"
#include "VectorFunctionDecimator.hpp"
#include <cxxtest/TestSuite.h>
//#include <iostream>

#include <vector>

class TestDecimation : public CxxTest::TestSuite
{
public:
    void TestBase1D() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1.0e-5);
        Decimator<1> base_decimator;
        base_decimator.Initialise(&mesh);
        TrianglesMeshWriter<1,1> mesh_writer("", "LineNoDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
        
        base_decimator.SetThreshold(0.5);
        base_decimator.Decimate();
        //base_decimator.Interrogate();
        
        TrianglesMeshWriter<1,1> mesh_writer1("", "LinePartDecimation");
        mesh_writer1.WriteFilesUsingMesh(mesh);
        //base_decimator.Interrogate();
        
        base_decimator.SetThreshold(2.0);
        base_decimator.Decimate();
        TrianglesMeshWriter<1,1> mesh_writer2("", "LineFullDecimation");
        mesh_writer2.WriteFilesUsingMesh(mesh);
        
    }
    void TestBase2DOnDisk()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        Decimator<2> base_decimator;
        base_decimator.Initialise(&mesh);
        
        TrianglesMeshWriter<2,2> mesh_writer("", "DiskNoDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
        base_decimator.SetThreshold(0.5);
        base_decimator.Decimate();
        
        //exit(1);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 113U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        TrianglesMeshWriter<2,2> mesh_writer1("", "DiskPartDecimation");
        mesh_writer1.WriteFilesUsingMesh(mesh);
        
        base_decimator.SetThreshold(INFINITY);
        base_decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 100U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        TrianglesMeshWriter<2,2> mesh_writer2("", "DiskFullDecimation");
        mesh_writer2.WriteFilesUsingMesh(mesh);
    }
    void TestBase2DOnSquare()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 441U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 800);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        Decimator<2> base_decimator;
        base_decimator.Initialise(&mesh);
        
        //base_decimator.Interrogate();
        TrianglesMeshWriter<2,2> mesh_writer("", "SquareNoDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
        
        base_decimator.SetThreshold(0.0005);
        base_decimator.Decimate();
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 80U);
        TS_ASSERT_LESS_THAN_EQUALS(73U, mesh.GetNumNodes());
        TrianglesMeshWriter<2,2> mesh_writer1("", "SquarePartDecimation");
        mesh_writer1.WriteFilesUsingMesh(mesh);
        
        
        //    base_decimator.Interrogate();
        base_decimator.SetThreshold(INFINITY);
        base_decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        TrianglesMeshWriter<2,2> mesh_writer2("", "SquareFullDecimation");
        mesh_writer2.WriteFilesUsingMesh(mesh);
        //     base_decimator.Interrogate();
    } 
    
    void TestBase2DOnSquareWithAutoMeshGeneration()
    {
       
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(20,20);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 400, 1.0e-5);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 441U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 800);
        Decimator<2> base_decimator;
        base_decimator.Initialise(&mesh);
        
        
        
        
        base_decimator.SetThreshold(6.000);
        base_decimator.DecimateAnimate("StructuredSquare");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 221U);
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 221U);
        TS_ASSERT_LESS_THAN_EQUALS(221U, mesh.GetNumNodes());
        
        
        base_decimator.SetThreshold(INFINITY);
        base_decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 400, 1.0e-5);
       
    }
    void TestBase3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1.0e-5);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 375U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1626);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 390);
        TrianglesMeshWriter<3,3> mesh_writer("", "CubeNoDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
        
        Decimator<3> base_decimator;
        base_decimator.Initialise(&mesh);
        base_decimator.SetThreshold(0.08);
        base_decimator.Decimate();
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 195U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 923);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 154);
        TrianglesMeshWriter<3,3> mesh_writer1("", "CubePartDecimation");
        mesh_writer1.WriteFilesUsingMesh(mesh);
        
        base_decimator.SetThreshold(INFINITY);
        base_decimator.Decimate();
        //  base_decimator.Interrogate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 8U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 12);
        TrianglesMeshWriter<3,3> mesh_writer2("", "CubeFullDecimation");
        mesh_writer2.WriteFilesUsingMesh(mesh);
    }
    void TestSequence2D()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        SequenceDecimator<2> decimator;
        decimator.Initialise(&mesh);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 441U);
        decimator.SetThreshold(102);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 341U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        decimator.SetThreshold(202);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 241U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        TrianglesMeshWriter<2,2> mesh_writer1("", "SquareSeqDecimation");
        mesh_writer1.WriteFilesUsingMesh(mesh);
        decimator.SetThreshold(302);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 141U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        decimator.SetThreshold(402);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 41U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        decimator.SetThreshold(440);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        
        decimator.SetThreshold(INFINITY);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
    }
    void TestSequence2DWithReIndex()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        SequenceDecimator<2> decimator;
        decimator.Initialise(&mesh);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 441U);
        decimator.SetThreshold(102);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 341U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        decimator.SetThreshold(202);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 241U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        
        mesh.ReIndex();
        TrianglesMeshWriter<2,2> mesh_writer1("", "SquareSeqDecimation");
        mesh_writer1.WriteFilesUsingMesh(mesh);

        decimator.Rescore();
        decimator.SetThreshold(102);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 141U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        decimator.SetThreshold(202);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 41U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        decimator.SetThreshold(240);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        
        decimator.SetThreshold(INFINITY);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
    }
    void TestSequence2DAnimate()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        SequenceDecimator<2> decimator;
        decimator.Initialise(&mesh);
        
        decimator.DecimateAnimate("SequentialAnimation");
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
    }
    void TestDiskQuality2D()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        QualityDecimator<2> decimator;
        decimator.Initialise(&mesh);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543U);
        decimator.SetThreshold(0.99);
        decimator.DecimateAnimate("DiskQuality");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 506U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        
        decimator.SetThreshold(INFINITY);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 506U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
    }   
    void TestQuality2D()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/SquareSeqDecimation");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        QualityDecimator<2> decimator;
        decimator.Initialise(&mesh);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 241U);
        decimator.SetThreshold(0.7);
        decimator.DecimateAnimate("Quality");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 187U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        
        decimator.SetThreshold(1.0);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 160U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        decimator.SetThreshold(INFINITY);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 160U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
    }   
    void TestMinimumElement2D()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        MinimumElementDecimator<2> decimator;
        decimator.Initialise(&mesh);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543U);
        //decimator.Interrogate();
        decimator.SetThreshold(1e-2);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 196U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        TrianglesMeshWriter<2,2> mesh_writer("", "DiskMimumumPartDecimation");
        mesh_writer.WriteFilesUsingMesh(mesh);
        
        decimator.SetThreshold(INFINITY);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 101U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        TrianglesMeshWriter<2,2> mesh_writer2("", "DiskMimumumFullDecimation");
        mesh_writer2.WriteFilesUsingMesh(mesh);
        
    }
    void TestMinimumElement2DAnimate()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        MinimumElementDecimator<2> decimator;
        decimator.Initialise(&mesh);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543U);
        decimator.DecimateAnimate("DiskMinimumAnimation");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 101U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
    }

    void Test1DLinearFunctionAnimate()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(),1.0,1e-6);
        
        //Make a payload
        std::vector<double> values(101);
        for (int i=0;i<101;i++)
        {
            values[i]=i;
        }
        LinearFunctionDecimator<1> decimator;
        decimator.Initialise(&mesh, values);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 101U);
        decimator.SetThreshold(0.0);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 101U);
        decimator.SetThreshold(2e-6);
        decimator.DecimateAnimate("LinearAnimation");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 2U);
        
        
    }
    void Test1DLinearFunctionSinAnimate()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(),1.0,1e-6);
        mesh.RescaleMeshFromBoundaryNode(Point<1>(2.0*M_PI),100);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(),6.2831853,1e-6);
        
        //Make a payload
        std::vector<double> values(101);
        for (int i=0;i<101;i++)
        {
            values[i]=sin(mesh.GetNode(i)->GetPoint()[0]);
        }
        LinearFunctionDecimator<1> decimator;
        decimator.Initialise(&mesh, values);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 101U);
        decimator.SetThreshold(0.0);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 101U);
        decimator.SetThreshold(0.01);
        decimator.DecimateAnimate("LinearSinAnimation");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 20U);
    }
    
    void Test2DLinearFunctionAnimate()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        //mesh.Translate(1.0,1.0);
   
        //Make a payload
        //This is binary state - a square in the middle has value 1 (like tumour)
        //                     - the surrounding disk has value -1 (like normal tissue)
        unsigned num_nodes=mesh.GetNumNodes();
        std::vector<double> values(num_nodes);
        for (unsigned i=0;i<num_nodes;i++)
        {
            double x=mesh.GetNode(i)->GetPoint()[0];
            double y=mesh.GetNode(i)->GetPoint()[1];
            if (fabs(x)<0.5 && fabs(y)<0.5)
            {
                values[i]=1.0;
            }
            else 
            {
                values[i]=-1.0;
            }
        }
        LinearFunctionDecimator<2> decimator;
        decimator.Initialise(&mesh, values);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);
        decimator.SetThreshold(0.001);
        decimator.DecimateAnimate("DiskLinearAnimation");
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 197U);
        TS_ASSERT_LESS_THAN_EQUALS(195U, mesh.GetNumNodes());
    }
    
    void Test2DVectorFunctionAnimate()
    {
        
        
        int width=40;
        int height=40;
        //int width=11;
        //int height=11;
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(width,height);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), width*height, 1.0e-5);
        mesh.Scale(2.0/width,2.0/height); //Now in [0,2]x[0,2]
        mesh.Translate(-1.0, -1.0); //Now in [-1,1]x[-1,1];
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 4.0, 1.0e-5);
        
        unsigned num_nodes=(width+1)*(height+1);
        std::vector<c_vector<double, 2> > values(num_nodes);
        double sqrt_two=sqrt(2.0);
        for (unsigned i=0;i<num_nodes;i++)
        {
            if (mesh.GetNode(i)->GetPoint()[0] >= 0.0)
            {
                values[i](0) = sqrt_two;
            }
            else 
            {
                values[i](0) = -sqrt_two;
            }
            if (mesh.GetNode(i)->GetPoint()[1] >= 0.0)
            {
                values[i](1) = sqrt_two;
            }
            else 
            {
                values[i](1) = -sqrt_two;
            }
        }
        
      
        VectorFunctionDecimator<2> decimator;
        decimator.Initialise(&mesh, values);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);
        decimator.SetThreshold(1e-6);
         decimator.DecimateAnimate("SquareVectorAnimation",10);
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 16U);
        
    }
    
};
#endif //_TESTDECIMATION_HPP_
