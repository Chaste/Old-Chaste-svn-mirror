#ifndef TESTCRYPTVORONOIDATAWRITER_HPP_
#define TESTCRYPTVORONOIDATAWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "FixedCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CryptVoronoiDataWriter.hpp"
#include "CellsGenerator.hpp"

class TestCryptVoronoiDataWriter : public CxxTest::TestSuite
{    
public:
    void TestDataWriter()
    {
        // set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        // create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 2, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();      
        
        // create the crypt
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, *p_mesh);
        
        Crypt<2> crypt(*p_mesh,cells);
        crypt.SetGhostNodes(ghost_node_indices);
        crypt.CreateVoronoiTessellation();
        
        // put this in brackets just so the writer does out scope, so its destructor
        // gets called and the file gets closed
        {        
            CryptVoronoiDataWriter<2> writer(crypt,"TestCryptVoronoiDataWriter","Simple.dat");
            writer.WriteData();
            p_simulation_time->IncrementTimeOneStep();
            writer.WriteData();
        }

        // work out where the previous test wrote its files
        OutputFileHandler handler("TestCryptVoronoiDataWriter",false);
        std::string results_file = handler.GetTestOutputDirectory() + "Simple.dat";
        TS_ASSERT_EQUALS(system(("ndiff -relerr 1e-5 " + results_file + " cancer/test/data/TestCryptVoronoiDataWriter/Simple.dat").c_str()), 0);
        
        SimulationTime::Destroy();
    }

    void TestDataWriterWithLoggedCell()
    {
        // set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        // create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 2, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();      
        
        // create the crypt
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, *p_mesh);
        
        Crypt<2> crypt(*p_mesh,cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        // flag a cell
        crypt.Begin()->SetLogged();
        
        crypt.CreateVoronoiTessellation();
        
        // put this in brackets just so the writer does out scope, so its destructor
        // gets called and the file gets closed
        {        
            CryptVoronoiDataWriter<2> writer(crypt,"TestCryptVoronoiDataWriter","OneCell.dat", true);
            writer.WriteData();
            p_simulation_time->IncrementTimeOneStep();
            writer.WriteData();
        }

        // work out where the previous test wrote its files
        OutputFileHandler handler("TestCryptVoronoiDataWriter",false);
        std::string results_file = handler.GetTestOutputDirectory() + "OneCell.dat";
        TS_ASSERT_EQUALS(system(("cmp " + results_file + " cancer/test/data/TestCryptVoronoiDataWriter/OneCell.dat").c_str()), 0);
        
        SimulationTime::Destroy();
    }
};


#endif /*TESTCRYPTVORONOIDATAWRITER_HPP_*/
