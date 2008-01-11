#ifndef TESTINGECELLSSIMULATION_HPP_
#define TESTINGECELLSSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "CellsGenerator.hpp"
#include "CryptSimulation2d.hpp"
#include "WntCellCycleModel.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "AbstractCellKiller.hpp"
#include "SloughingCellKiller.hpp"

class TestIngeCellsSimulation : public CxxTest::TestSuite
{
public:
    // Test to check that none of the protein concentrations become -ve for different Wnt stimuli
    // See ticket 629 and bodge in IngeWntSwatCellCycleModel.cpp in initial conditions
    // Note this is fairly random on which Wnt stimuli make protein conc 9 and 17 go -ve.
    // !!Don't change the crypt height on this test unless you know that the bodge is being used.!! 
    
    void TestIngeBetaCatVis() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
                
        //double end_of_simulation = 150.0; // hours
        double time_of_each_run = 0.01; // for each run
        
        unsigned cells_across = 4;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, INGE_WNT_SWAT_HYPOTHESIS_ONE, true);
        
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i].SetBirthTime(-1.1); // Just to make the test run a bit quicker.
        }
        
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        WntGradient::Instance()->SetType(OFFSET_LINEAR);
        WntGradient::Instance()->SetTissue(crypt);
        
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("IngeCellsNiceCryptSim_long");
        
        // Set simulation to output cell types
        simulator.SetOutputCellTypes(true);
                
        // Set length of simulation here
        simulator.SetEndTime(time_of_each_run);
        
        simulator.SetMaxCells(1000);
        simulator.SetMaxElements(2000);

        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&simulator.rGetTissue(),0.01);
        simulator.AddCellKiller(p_cell_killer);
        
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        simulator.Save();
                
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
};




#endif /*TESTINGECELLSSIMULATION_HPP_*/
