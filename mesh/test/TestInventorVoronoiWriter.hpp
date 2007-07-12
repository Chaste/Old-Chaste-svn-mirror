#ifndef TESTINVENTORVORONOIWRITER_HPP_
#define TESTINVENTORVORONOIWRITER_HPP_


#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"
#include "VoronoiTessellation.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Exception.hpp"
#include "InventorVoronoiWriter.hpp"
#include "TrianglesMeshReader.cpp"

#include <cmath>
#include <vector>

class TestInventorVoronoiWriter : public CxxTest::TestSuite
{
public:
    void TestWriteTet() throw (Exception)
    {
        
        // Create conforming tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;
        
        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,0.5,0.5));
        
        ConformingTetrahedralMesh<3,3> mesh(nodes);
        
        // Create Voronoi Tesselation
        VoronoiTessellation tessellation(mesh);
        
        InventorVoronoiWriter inventor_writer("InventorWriter", "SimpleTet");
        inventor_writer.Write(tessellation);
        
        // then compare against known good file
        OutputFileHandler handler("InventorWriter",false);
        std::string results_dir = handler.GetTestOutputDirectory();
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/SimpleTet.iv mesh/test/data/InventorWriter/SimpleTet.iv").c_str()), 0);
        
    }
    
    void TestWriteComplexCube() throw (Exception)
    {
        
        // Create conforming tetrahedral mesh which is Delauny
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_136_elements");
        
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        TS_ASSERT(mesh.CheckVoronoi());
        
        // Create Voronoi Tesselation
        VoronoiTessellation tessellation(mesh);
        
        InventorVoronoiWriter inventor_writer("InventorWriter", "Complex", false);
        inventor_writer.Write(tessellation);
        
        // then compare against known good file
        OutputFileHandler handler("InventorWriter",false);
        std::string results_dir = handler.GetTestOutputDirectory();
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/Complex.iv mesh/test/data/InventorWriter/Complex.iv").c_str()), 0);
        
    }
};
#endif /*TESTINVENTORVORONOIWRITER_HPP_*/
