#ifndef TESTMAKENICECRYPTSIMS_HPP_
#define TESTMAKENICECRYPTSIMS_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "TissueCell.hpp"
#include "WntCellCycleModel.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "AbstractCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "CellsGenerator.hpp"

class TestMakeNiceCryptSims : public CxxTest::TestSuite
{
public:
void TestNiceCryptSimulation() throw (Exception)
    {        
        CancerParameters* p_params = CancerParameters::Instance();
        std::string output_directory = "NiceCryptSim";
        
        double end_of_simulation = 500.0; // hours
        double time_of_each_run = 50.0; // for each run
        
        unsigned cells_across = 13;
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
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, STOCHASTIC_WNT, true);
              
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
                
        WntGradient::Instance()->SetType(OFFSET_LINEAR);
        WntGradient::Instance()->SetCrypt(crypt);
        
        TissueSimulation<2> simulator(crypt);
        simulator.SetOutputDirectory(output_directory);
        
        // Set simulation to output cell types
        simulator.SetOutputCellTypes(true);
                
        // Set length of simulation here
        simulator.SetEndTime(time_of_each_run);
        
        simulator.SetMaxCells(1000);
        simulator.SetMaxElements(2000);

        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&simulator.rGetTissue(),0.01);
        simulator.AddCellKiller(p_cell_killer);
        
        // UNUSUAL SET UP HERE /////////////////////////////////////
        
        p_params->SetDampingConstantNormal(1.0);    // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());
        
        p_params->SetSpringStiffness(30.0); //normally 15.0;
        // 0.3/30 = 0.01 (i.e. Meineke's values)
        
        simulator.UseNonFlatBottomSurface();
        
        // END OF UNUSUAL SET UP! //////////////////////////////////
        
        simulator.Solve();
        
        std::vector<unsigned> label_these = Label();
        // set a cell to be labelled (probably a stemish cell)
        for (unsigned i=0; i<label_these.size() ; i++)
        {
            simulator.rGetTissue().rGetCellAtNodeIndex(label_these[i]).SetMutationState(LABELLED);
        }
        
        simulator.Save();
        
        for (double t=time_of_each_run; t<end_of_simulation+0.5; t += time_of_each_run)
        {
            TissueSimulation<2>* p_simulator = TissueSimulation<2>::Load("NiceCryptSim",t);
            p_simulator->SetEndTime(t+time_of_each_run);
            p_simulator->Solve();
            p_simulator->Save();
            delete p_simulator;
        }
                
        delete p_cell_killer;
        delete p_params;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
    
std::vector<unsigned> Label()
{
    std::vector<unsigned> label_these;
//    label_these.push_back(37);
//    label_these.push_back(309);
    return label_these;   
}

};

#endif /*TESTMAKENICECRYPTSIMS_HPP_*/
