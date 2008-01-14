#ifndef TESTCRYPTVORONOIDATAWRITER_HPP_
#define TESTCRYPTVORONOIDATAWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "FixedCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TissueVoronoiDataWriter.hpp"
#include "CellsGenerator.hpp"

/**
 * Note that all these tests call setUp() and tearDown() before running,
 * so if you copy them into a new test suite be sure to copy these methods
 * too.
 */
class TestTissueVoronoiDataWriter : public CxxTest::TestSuite
{    
private:

    void setUp()
    {
        // Initialise singleton classes
        SimulationTime::Instance()->SetStartTime(0.0);
    }
    void tearDown()
    {
        // Clear up singleton classes
        SimulationTime::Destroy();
    }
    
public:

    void TestDataWriter()
    {
        // Set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 2, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();      
        
        // Create the crypt
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, *p_mesh);
        
        Tissue<2> tissue(*p_mesh,cells);
        tissue.SetGhostNodes(ghost_node_indices);
        tissue.CreateVoronoiTessellation();
        
        // Put this in brackets just so the writer does out scope, so its destructor
        // gets called and the file gets closed
        {        
            TissueVoronoiDataWriter<2> writer(tissue,"TestTissueVoronoiDataWriter","Simple.dat");
            writer.WriteData();
            p_simulation_time->IncrementTimeOneStep();
            writer.WriteData();
        }

        // Work out where the previous test wrote its files
        OutputFileHandler handler("TestTissueVoronoiDataWriter",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "Simple.dat";
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/TestTissueVoronoiDataWriter/Simple.dat").c_str()), 0);
    }

    void TestDataWriterWithLoggedCell()
    {
        // Set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 2, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();      
        
        // Create the crypt
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, *p_mesh);
        
        Tissue<2> tissue(*p_mesh,cells);
        tissue.SetGhostNodes(ghost_node_indices);
        
        // Flag a cell
        tissue.Begin()->SetLogged();
        
        tissue.CreateVoronoiTessellation();
        
        // Put this in brackets just so the writer goes out of scope, 
        // so its destructor gets called and the file gets closed.
        {        
            TissueVoronoiDataWriter<2> writer(tissue,"TestTissueVoronoiDataWriter","OneCell.dat", true);
            writer.WriteData();
            p_simulation_time->IncrementTimeOneStep();
            writer.WriteData();
        }

        // Work out where the previous test wrote its files
        OutputFileHandler handler("TestTissueVoronoiDataWriter",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "OneCell.dat";
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/TestTissueVoronoiDataWriter/OneCell.dat").c_str()), 0);
    }
    
    void TestWriteTissueAreas()
    {
        // Set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create a simple mesh
        HoneycombMeshGenerator generator(4, 4, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();   

        // Create the crypt
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, *p_mesh);

        // Set one of the non-boundary cells to be necrotic
        cells[6].SetCellType(NECROTIC);
        
        Tissue<2> tissue(*p_mesh,cells);
        tissue.SetGhostNodes(ghost_node_indices);
        tissue.CreateVoronoiTessellation();

        // Put this in brackets just so the writer does out scope, so its destructor
        // gets called and the file gets closed
        {        
            TissueVoronoiDataWriter<2> writer(tissue,"TestTissueVoronoiDataWriter","Areas.dat");
            writer.WriteTissueAreas();
            p_simulation_time->IncrementTimeOneStep();
            writer.WriteTissueAreas();
        }

        // Work out where the previous test wrote its files
        OutputFileHandler handler("TestTissueVoronoiDataWriter",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "Areas.dat";
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/TestTissueVoronoiDataWriter/Areas.dat").c_str()), 0);

        // Expect 3.4641016 for total area and 0.8660255 for necrotic
        // because each cell is 0.8660255 and there is one nectoric and four non-boundary cells in total.
    }
};


#endif /*TESTTISSUEVORONOIDATAWRITER_HPP_*/
