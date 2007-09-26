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

class TestCryptVoronoiDataWriter : public CxxTest::TestSuite
{    
private: 
///\todo: this is identical to SetUpCells in crypt - refactor
    template<unsigned DIM>
    std::vector<MeinekeCryptCell> SetUpCells(ConformingTetrahedralMesh<DIM,DIM>* pMesh)
    {
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<pMesh->GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        return cells;
    }
    
public:
    void TestDataWriter()
    {
        // set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 2, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();      
        
        // create the crypt
        std::vector<MeinekeCryptCell> cells = SetUpCells<2>(p_mesh);
        Crypt<2> crypt(*p_mesh,cells);
        crypt.SetGhostNodes(ghost_node_indices);
        crypt.CreateVoronoiTessellation();
        
        // put this in brackets just so the writer does out scope, so it's destructor
        // gets called and the file in closed
        {        
            CryptVoronoiDataWriter writer(crypt,"TestCryptVoronoiDataWriter","Simple.dat");
            writer.WriteData();
        }

        // work out where the previous test wrote its files
        OutputFileHandler handler("TestCryptVoronoiDataWriter",false);
        std::string results_file = handler.GetTestOutputDirectory() + "Simple.dat";
        TS_ASSERT_EQUALS(system(("cmp " + results_file + " cancer/test/data/TestCryptVoronoiDataWriter/Simple.dat").c_str()), 0);
        
        SimulationTime::Destroy();
    }
};


#endif /*TESTCRYPTVORONOIDATAWRITER_HPP_*/
