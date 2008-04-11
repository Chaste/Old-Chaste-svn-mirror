#ifndef TESTGENERATESTEADYSTATECRYPT_HPP_
#define TESTGENERATESTEADYSTATECRYPT_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <vector>

#include "CryptSimulation2d.hpp"
#include "WntCellCycleModel.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "TrianglesMeshReader.hpp"
#include "OutputFileHandler.hpp"


class TestGenerateSteadyStateCrypt : public CxxTest::TestSuite
{
public:

    /*
     * This test can be used to generate steady state crypts for use
     * as the starting points of other simulations.
     * 
     * You need to specify :
     * the kind of cell cycle model to use on line 64,
     * WntConcentration on line 69,
     * change any cancer parameters around line 90,
     * and give the simulator options around line 95.
     */
    void TestGenerateSteadyStateCryptArchives() throw (Exception)
    {        
        CancerParameters* p_params = CancerParameters::Instance();
        std::string output_directory = "SteadyStateCrypt";
        
        double end_of_simulation = 150.0; // hours
        double time_of_each_run = 50.0; // for each run
        
        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, STOCHASTIC_WNT, true);
              
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        WntConcentration::Instance()->SetType(LINEAR);
        CancerParameters::Instance()->SetTopOfLinearWntConcentration(1.0/3.0);
        WntConcentration::Instance()->SetTissue(crypt);
        
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory(output_directory);
        
        // Set simulation to output cell types
        simulator.SetOutputCellMutationStates(true);
                
        // Set length of simulation here
        simulator.SetEndTime(time_of_each_run);

        SloughingCellKiller cell_killer(&simulator.rGetTissue(),0.01);
        simulator.AddCellKiller(&cell_killer);
        
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
        
        for (double t=time_of_each_run; t<end_of_simulation+0.5; t += time_of_each_run)
        {
            CryptSimulation2d* p_simulator = CryptSimulation2d::Load("SteadyStateCrypt",t);
            p_simulator->SetEndTime(t+time_of_each_run);
            p_simulator->Solve();
            p_simulator->Save();
            delete p_simulator;
        }
                
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntConcentration::Destroy();
    }
    
};

#endif /*TESTGENERATESTEADYSTATECRYPT_HPP_*/
