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

class TestRepresentativeSimulation : public CxxTest::TestSuite
{
public:
void TestIngeBetaCatVis() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        
        //double end_of_simulation = 150.0; // hours
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
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, INGE_WNT_SWAT_HYPOTHESIS_ONE, true);
              
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
        
        // UNUSUAL SET UP HERE /////////////////////////////////////
        
        p_params->SetDampingConstantNormal(1.0);    // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());
        
        p_params->SetSpringStiffness(30.0); //normally 15.0;
        // 0.3/30 = 0.01 (i.e. Meineke's values)
        
        simulator.UseJiggledBottomCells();
        
        // END OF UNUSUAL SET UP! //////////////////////////////////
        
                
        simulator.Solve();
        simulator.Save();
//        double end_of_simulation = 350.0; // hours
//        double time_of_each_run = 50.0; // for each run
//        double start_time = 200.0;
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(50.0);
//        
//        
//        for (double t=time_of_each_run; t<end_of_simulation+0.5; t += time_of_each_run)
//        {
//            CryptSimulation2d* p_simulator = CryptSimulation2d::Load("IngeCellsNiceCryptSim_long",t);
//            p_simulator->SetEndTime(t+time_of_each_run);
//            p_simulator->Solve();
//            p_simulator->Save();
//            delete p_simulator;
//        }
                
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
std::vector<unsigned> Label()
{
    std::vector<unsigned> label_these;
//    label_these.push_back(442);
//    label_these.push_back(290);
//    label_these.push_back(417);
//    label_these.push_back(260);
    label_these.push_back(505);
    label_these.push_back(206);
    label_these.push_back(40);
    return label_these;   
}
};




#endif /*TESTINGECELLSSIMULATION_HPP_*/
