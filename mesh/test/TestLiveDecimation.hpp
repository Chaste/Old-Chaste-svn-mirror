#ifndef _TESTLIVEDECIMATION_HPP_
#define _TESTLIVEDECIMATION_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include "LinearFunctionDecimator.hpp"
#include "VectorFunctionDecimator.hpp"
#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>
#include <fstream>

class TestLiveDecimation : public CxxTest::TestSuite
{
public:
    void TestCancerExample() throw (Exception)
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
            if (data == '\x00'){
                values[i]=-1.0;
            } else {
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
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumElements(), 810);
        TS_ASSERT_LESS_THAN_EQUALS(786, mesh.GetNumElements());
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumBoundaryElements(), 12);
        TS_ASSERT_LESS_THAN_EQUALS(11, mesh.GetNumBoundaryElements());
        
    }
    void TestJoeExample() throw (Exception)
    {
        std::ifstream input_ppm("mesh/test/data/joe.pgm", std::ios::binary);
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
            if (data == '\x00'){
                values[i]=-1.0;
            } else {
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
        decimator.DecimateAnimate("JoeAnimation", 15);
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 301U);
        TS_ASSERT_LESS_THAN_EQUALS(297U, mesh.GetNumNodes());
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumElements(), 588);
        TS_ASSERT_LESS_THAN_EQUALS(581, mesh.GetNumElements());
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumBoundaryElements(), 14);
        TS_ASSERT_LESS_THAN_EQUALS(10, mesh.GetNumBoundaryElements());
    }
    void TestEndExample() throw (Exception)
    {
        std::ifstream input_ppm("mesh/test/data/end.pgm", std::ios::binary);
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
            if (data == '\x00'){
                values[i]=-1.0;
            } else {
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
        decimator.DecimateAnimate("EndAnimation", 100);
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumNodes(), 346U);
        TS_ASSERT_LESS_THAN_EQUALS(346U, mesh.GetNumNodes());
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumElements(), 679);
        TS_ASSERT_LESS_THAN_EQUALS(679, mesh.GetNumElements());
        TS_ASSERT_LESS_THAN_EQUALS(mesh.GetNumBoundaryElements(), 14);
        TS_ASSERT_LESS_THAN_EQUALS(10, mesh.GetNumBoundaryElements());
    }
    
    void TestTahirExample() throw (Exception)
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
#endif //_TESTLIVEDECIMATION_HPP_
