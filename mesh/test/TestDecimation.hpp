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
#include "RandomDecimator.hpp"
#include <cxxtest/TestSuite.h>
//#include <iostream>

#include <vector>

template <unsigned SPACE_DIM>
class PositionDecimator : public Decimator<SPACE_DIM>
{
protected:

    void CalculateLocalMeasure(NodeInfo<SPACE_DIM> *pNodeInfo, bool before)
    {
        if (before)
        {
            double measure=pNodeInfo->pGetNode()->rGetLocation()[0];//x-value
            this->mMeasureBefore=measure;
        }
    }
    double CalculateScore()
    {
        return this->mMeasureBefore;
    }
      

};

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
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 11U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1.0e-5);
        //TrianglesMeshWriter<1,1> mesh_writer("", "LineNoDecimation");
        //mesh_writer.WriteFilesUsingMesh(mesh);
        
        base_decimator.SetThreshold(0.5);
        base_decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1.0e-5);
        //base_decimator.Interrogate();
        
        //TrianglesMeshWriter<1,1> mesh_writer1("", "LinePartDecimation");
        //mesh_writer1.WriteFilesUsingMesh(mesh);
        //base_decimator.Interrogate();
        
        base_decimator.SetThreshold(2.0);
        base_decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 2U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1.0e-5);
        //TrianglesMeshWriter<1,1> mesh_writer2("", "LineFullDecimation");
        //mesh_writer2.WriteFilesUsingMesh(mesh);
        
    }
    void TestBase2DOnDisk()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        Decimator<2> base_decimator;
        base_decimator.Initialise(&mesh);
        
        //TrianglesMeshWriter<2,2> mesh_writer("", "DiskNoDecimation");
        //mesh_writer.WriteFilesUsingMesh(mesh);
        base_decimator.SetThreshold(0.5);
        base_decimator.Decimate();
        
        
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 113U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        //TrianglesMeshWriter<2,2> mesh_writer1("", "DiskPartDecimation");
        //mesh_writer1.WriteFilesUsingMesh(mesh);
        
        base_decimator.SetThreshold(INFINITY);
        base_decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 100U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        //TrianglesMeshWriter<2,2> mesh_writer2("", "DiskFullDecimation");
        //mesh_writer2.WriteFilesUsingMesh(mesh);
    }
    
    void TestBase2DOnDiskWithVolumeLeak()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.13953, 1.0e-5);
        Decimator<2> decimator;
        decimator.Initialise(&mesh);
        
        
        
        TS_ASSERT_DELTA(decimator.GetVolumeLeakage(), 1e-5, 1.0e-10);
        decimator.SetVolumeLeakage(1e-2);
        
        decimator.SetThreshold(0.5);
        decimator.Decimate();
        
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 46U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.1133, 1.0e-4);
        
        
        QualityDecimator<2> quality_decimator;
        quality_decimator.Initialise(&mesh);
        quality_decimator.SetThreshold(0.3);
        quality_decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 41U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 3.1133, 1.0e-4);
        
        //TrianglesMeshWriter<2,2> mesh_writer1("", "DiskLeakage");
        //mesh_writer1.WriteFilesUsingMesh(mesh);
        
    }
    
    void TestBase2DOnSquare()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 441U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 800U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        Decimator<2> base_decimator;
        base_decimator.Initialise(&mesh);
        
        //base_decimator.Interrogate();
        //TrianglesMeshWriter<2,2> mesh_writer("", "SquareNoDecimation");
        //mesh_writer.WriteFilesUsingMesh(mesh);
        
        base_decimator.SetThreshold(0.0005);
        base_decimator.Decimate();
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 80U);
        TS_ASSERT_LESS_THAN_EQUALS(73U, mesh.GetNumNodes());
        //TrianglesMeshWriter<2,2> mesh_writer1("", "SquarePartDecimation");
        //mesh_writer1.WriteFilesUsingMesh(mesh);
        
        
        //    base_decimator.Interrogate();
        base_decimator.SetThreshold(INFINITY);
        base_decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        //TrianglesMeshWriter<2,2> mesh_writer2("", "SquareFullDecimation");
        //mesh_writer2.WriteFilesUsingMesh(mesh);
        //     base_decimator.Interrogate();
    }
    
    void TestBase2DOnSquareWithStopping()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 441U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 800U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        FixedNodeDecimator<2> forty_node_decimator(40);
        
        forty_node_decimator.Initialise(&mesh);
        
         
        forty_node_decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 40U);
      
    }
    
       void TestPositionOnSquare()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 441U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 800U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        PositionDecimator<2> decimator;
        
        decimator.Initialise(&mesh);
        TS_ASSERT_DELTA(decimator.GetVolumeLeakage(), 1e-5, 1.0e-10);
        decimator.SetVolumeLeakage(100);
        decimator.SetThreshold(0.05);
         
        //decimator.Decimate();
        decimator.DecimateAnimate("PositionDecimator");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 232U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 420U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.0052, 1.0e-4);
      
    }
    
    
    void TestRandom2DOnSquareWithAutoMeshGeneration()
    {
    
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(20,20);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 400, 1.0e-5);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 441U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 800U);
        RandomDecimator<2> decimator;
        RandomNumberGenerator::Instance();
        decimator.Initialise(&mesh);
        
        
        //decimator.Interrogate();
        
        decimator.SetThreshold(0.2);
        decimator.DecimateAnimate("RandomSquare");
        //TS_ASSERT_EQUALS(mesh.GetNumNodes(), 254U);
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 263U);
        TS_ASSERT_LESS_THAN_EQUALS(254U, mesh.GetNumNodes());
        
        
        ConformingTetrahedralMesh<2,2> mesh2;
        mesh2.ConstructRectangularMesh(20,20);
        RandomDecimator<2> decimator2;
        decimator2.Initialise(&mesh2);
        
        
        //decimator.Interrogate();
        
        decimator2.SetThreshold(0.2);
        decimator2.Decimate();
        //TS_ASSERT_EQUALS(mesh2.GetNumNodes(), 254U);
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 263U);
        TS_ASSERT_LESS_THAN_EQUALS(254U, mesh.GetNumNodes());
        
        //Note that each re-score gives a 50% chance of a node falling below 0.5
        decimator2.SetThreshold(0.8);
        decimator2.Decimate();
        TS_ASSERT_LESS_THAN_EQUALS(mesh2.GetNumNodes(), 6U);
        TS_ASSERT_LESS_THAN_EQUALS(4U, mesh2.GetNumNodes());
        TS_ASSERT_DELTA(mesh2.CalculateMeshVolume(), 400, 1.0e-5);
        RandomNumberGenerator::Destroy();
        
    }
    void TestBase3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1.0e-5);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 375U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1626U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 390U);
        //TrianglesMeshWriter<3,3> mesh_writer("", "CubeNoDecimation");
        //mesh_writer.WriteFilesUsingMesh(mesh);
        
        Decimator<3> base_decimator;
        base_decimator.Initialise(&mesh);
        base_decimator.SetThreshold(0.08);
        base_decimator.Decimate();
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 189U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 906U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 144U);
        //TrianglesMeshWriter<3,3> mesh_writer1("", "CubePartDecimation");
        //mesh_writer1.WriteFilesUsingMesh(mesh);
        
        base_decimator.SetThreshold(INFINITY);
        base_decimator.Decimate();
        //  base_decimator.Interrogate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 8U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 12U);
        //TrianglesMeshWriter<3,3> mesh_writer2("", "CubeFullDecimation");
        //mesh_writer2.WriteFilesUsingMesh(mesh);
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
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 179U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        
        decimator.SetThreshold(1.0);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 156U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
        decimator.SetThreshold(INFINITY);
        decimator.Decimate();
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 156U);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01, 1.0e-5);
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
    
    void Test1DLinearFunction()
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
        decimator.Decimate();
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
    
    
    void xTestCancerExample() throw (Exception)
    {
        std::ifstream input_ppm("mesh/test/data/filled_colon_downsized.pgm", std::ios::binary);
        std::string magic;
        input_ppm >> magic;
        
        TS_ASSERT_EQUALS(magic.compare("P5"), 0);
        
        unsigned x, y, levels;
        input_ppm >> x >> y >> levels;
        TS_ASSERT_EQUALS(levels, 255U);
        
        char data;
        input_ppm.read(&data,1);
        TS_ASSERT_EQUALS(data,'\n');
        
        //Make a payload
        //This is binary state
        unsigned num_nodes=x*y;
        std::vector<double> values(num_nodes);
        for (unsigned i=0;i<num_nodes;i++)
        {
            input_ppm.read(&data,1);
            if (data == '\x00')
            {
                values[i]=-1.0;
            }
            else
            {
                TS_ASSERT_EQUALS(data, '\xff');
                values[i]=1.0;
            }
            
            
        }
        
        //Make a mesh
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(x-1,y-1);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), (x-1)*(y-1), 1.0e-5);
        
        
        
        LinearFunctionDecimator<2> decimator;
        decimator.Initialise(&mesh, values);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);
        decimator.SetThreshold(0.001);
        decimator.DecimateAnimate("ColonAnimation", 15);
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 412U);
        TS_ASSERT_LESS_THAN_EQUALS(400U, mesh.GetNumNodes());
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumElements(), 810U);
        TS_ASSERT_LESS_THAN_EQUALS(786U, mesh.GetNumElements());
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumBoundaryElements(), 12U);
        TS_ASSERT_LESS_THAN_EQUALS(11U, mesh.GetNumBoundaryElements());
        
    }
    
    void xTestTahirExample() throw (Exception)
    {
        std::ifstream input_fibres("mesh/test/data/fibres/Hist_Fibre_Vectors.txt");
        
        
        unsigned width=19;
        unsigned height=19;
        
        unsigned num_nodes=width*height;
        std::vector<c_vector<double, 2> > values(num_nodes);
        //Read x values
        for (unsigned i=0;i<num_nodes;i++)
        {
            double data;
            input_fibres>>data;
            values[i](0) = data;
        }
        //Read y values
        for (unsigned i=0;i<num_nodes;i++)
        {
            double data;
            input_fibres>>data;
            values[i](1) = data;
        }
        
        //Normalise directions
        for (unsigned i=0;i<num_nodes;i++)
        {
            double norm=norm_2(values[i]);
            
            TS_ASSERT_LESS_THAN(0.3, norm);
            values[i]=values[i]/norm;
        }
        for (unsigned i=0;i<num_nodes;i++)
        {
        
            TS_ASSERT_DELTA(norm_2(values[i]), 1.0, 1e-7);
        }
        
        //Make a mesh
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(width-1,height-1);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), (width-1)*(height-1), 1.0e-5);
        
        
        
        VectorFunctionDecimator<2> decimator;
        decimator.Initialise(&mesh, values);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);
        decimator.SetThreshold(2.0);
        decimator.DecimateAnimate("Tahir", 10);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 60U);
        //decimator.Interrogate();
        
    }
    
};
#endif //_TESTDECIMATION_HPP_
