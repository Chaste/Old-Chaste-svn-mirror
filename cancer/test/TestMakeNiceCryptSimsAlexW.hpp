#ifndef TESTMAKENICECRYPTSIMSALEXW_HPP_
#define TESTMAKENICECRYPTSIMSALEXW_HPP_

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

class TestMakeNiceCryptSimsAlexW : public CxxTest::TestSuite
{
public:
//void xTestRepresentativeSimulationForProfiling() throw (Exception)
//    {        
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//
//        std::string test_to_profile = "NiceCryptSim";
//        double t = 350;   // this is the folder and time that the stored results were archived (needed to know foldernames)
//        double run_for = 10; // run for 10 hours.
//        
//        // The archive needs to be copied from cancer/test/data/<test_to_profile>
//        // to the testoutput directory to continue running the simulation.     
//        OutputFileHandler any_old_handler("",false);
//        std::string test_output_directory = any_old_handler.GetTestOutputDirectory();
//        std::string test_data_directory = "cancer/test/data/" + test_to_profile +"/";
//        std::string command = "cp -R --remove-destination --reply=yes " + test_data_directory +" "+ test_output_directory +"/";     
//        int return_value = system(command.c_str());
//        TS_ASSERT_EQUALS(return_value, 0);
//        
//        TissueSimulation<2>* p_simulator = TissueSimulation<2>::Load(test_to_profile,t);
//        
//        std::vector<unsigned> label_these = Label();
//        // set a cell to be labelled (probably a stemish cell)
//        for (unsigned i=0; i<label_these.size() ; i++)
//        {
//            p_simulator->rGetCrypt().rGetCellAtNodeIndex(label_these[i]).SetMutationState(APC_TWO_HIT);
//        }
//        
//         // UNUSUAL SET UP HERE /////////////////////////////////////
//        CancerParameters* p_params = CancerParameters::Instance();
//        p_params->SetDampingConstantNormal(1.0);    // normally 1
//
//        // Do not give mutant cells any different movement properties to normal ones
//        p_params->SetDampingConstantMutant(10.0);
//        
//        p_params->SetSpringStiffness(30.0); //normally 15.0;
//        // 0.3/30 = 0.01 (i.e. Meineke's values)
//        
//        p_simulator->UseNonFlatBottomSurface();
//        
//        // END OF UNUSUAL SET UP! //////////////////////////////////
//        
//        p_simulator->SetEndTime(t+50); // start time + duration
//        p_simulator->Solve();
//        p_simulator->Save();
//        delete p_simulator;
//        
//        TissueSimulation<2>* p_simulator2 = TissueSimulation<2>::Load("NiceCryptSim",400);
//        p_simulator2->SetEndTime(450);
//        p_simulator2->Solve();
//        p_simulator2->Save();
//        delete p_simulator2;
//        
//        TissueSimulation<2>* p_simulator3 = TissueSimulation<2>::Load("NiceCryptSim",450);
//        p_simulator3->SetEndTime(475);
//        p_simulator3->Solve();
//        p_simulator3->Save();
//        delete p_simulator3;
//                
//        SimulationTime::Destroy();
//        RandomNumberGenerator::Destroy();
//    }
    
//void TestIngeBetaCatVis() throw (Exception)
//    {
//        CancerParameters *p_params = CancerParameters::Instance();
//        p_params->Reset();
//        // There is no limit on transit cells in Wnt simulation
//        p_params->SetMaxTransitGenerations(1000);
//        
//        
//        double time_of_each_run = 30.0; // for each run
//        
//        unsigned cells_across = 13;
//        unsigned cells_up = 25;
//        double crypt_width = 12.1;
//        unsigned thickness_of_ghost_layer = 3;
//        
//        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
//        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
//        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
//        
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        
//        // Set up cells
//        std::vector<TissueCell> cells;
//        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, INGE_WNT_SWAT_HYPOTHESIS_TWO, true);
//              
//        Tissue<2> crypt(*p_mesh, cells);
//        crypt.SetGhostNodes(ghost_node_indices);
//                
//        WntGradient::Instance()->SetType(LINEAR);
//        CancerParameters::Instance()->SetTopOfLinearWntGradient(1.0/3.0);
//        WntGradient::Instance()->SetTissue(crypt);
//        
//        CryptSimulation2d simulator(crypt);
//        simulator.SetOutputDirectory("IngeCellsNiceCryptSim_hyp2_long");
//        
//        // Set simulation to output cell types
//        simulator.SetOutputCellTypes(true);
//                
//        // Set length of simulation here
//        simulator.SetEndTime(time_of_each_run);
//
//        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&simulator.rGetTissue(),0.01);
//        simulator.AddCellKiller(p_cell_killer);
//        
//        // UNUSUAL SET UP HERE /////////////////////////////////////
//        
//        p_params->SetDampingConstantNormal(1.0);    // normally 1
//
//        // Do not give mutant cells any different movement properties to normal ones
//        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());
//        
//        p_params->SetSpringStiffness(30.0); //normally 15.0;
//        // 0.3/30 = 0.01 (i.e. Meineke's values)
//        
//        simulator.UseJiggledBottomCells();
//        
//        // END OF UNUSUAL SET UP! //////////////////////////////////
//        
//        std::cout<< "About to solve \n" << std::flush;
//        simulator.Solve();
//        simulator.Save();
//        double end_of_simulation = 350.0; // hours
//        
//        std::cout<< "Going into loop \n" << std::flush;
//        
//        for (double t=time_of_each_run; t<end_of_simulation+0.5; t += time_of_each_run)
//        {
//            std::cout<< "Results from time " << time_of_each_run << "\n" << std::flush;
//            CryptSimulation2d* p_simulator = CryptSimulation2d::Load("IngeCellsNiceCryptSim_hyp2_long",t);
//            p_simulator->SetEndTime(t+time_of_each_run);
//            p_simulator->Solve();
//            p_simulator->Save();
//            delete p_simulator;
//        }
//                
//        delete p_cell_killer;
//        SimulationTime::Destroy();
//        RandomNumberGenerator::Destroy();
//        WntGradient::Destroy();
//    }
void TestAreaDependentAndLengthDependent() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        
        double time_of_each_run = 3.0; // for each run
        
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
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, SIMPLE_WNT, true);
              
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
                
        WntGradient::Instance()->SetType(LINEAR);
        CancerParameters::Instance()->SetTopOfLinearWntGradient(1.0/3.0);
        WntGradient::Instance()->SetTissue(crypt);
        
        
        Meineke2001SpringSystem<2> meineke_spring_system(crypt);
        meineke_spring_system.SetAreaBasedViscosity(true);
        meineke_spring_system.SetEdgeBasedSpringConstant(true);
        
        CryptSimulation2d simulator(crypt, &meineke_spring_system, false, true);
        simulator.SetOutputDirectory("Noddy_WNT_Yes_Area_Yes_Length");
        
        // Set simulation to output cell types
        simulator.SetOutputCellTypes(true);
                
        // Set length of simulation here
        simulator.SetEndTime(time_of_each_run);

        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&simulator.rGetTissue(),0.01);
        simulator.AddCellKiller(p_cell_killer);
        
        // UNUSUAL SET UP HERE /////////////////////////////////////
        
        p_params->SetDampingConstantNormal(1.0);    // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());
        
        p_params->SetSpringStiffness(30.0); //normally 15.0;
        // 0.3/30 = 0.01 (i.e. Meineke's values)
        
        simulator.UseJiggledBottomCells();
        simulator.SetWriteVoronoiData(true, false); //Writes Area and perimeter
        
        // END OF UNUSUAL SET UP! //////////////////////////////////
        std::cout<< "About to solve \n" << std::flush;
        simulator.Solve();
//        simulator.Save();
//        double end_of_simulation = 400.0; // hours
//        
//        std::cout<< "Going into loop \n" << std::flush;
//        
//        for (double t=time_of_each_run; t<end_of_simulation+0.5; t += time_of_each_run)
//        {
//            std::cout<< "Results from time " << time_of_each_run << "\n" << std::flush;
//            CryptSimulation2d* p_simulator = CryptSimulation2d::Load("Noddy_WNT_Yes_Area_Yes_Length",t);
//            p_simulator->SetEndTime(t+time_of_each_run);
//            p_simulator->Solve();
//            p_simulator->Save();
//            delete p_simulator;
//        }
                
        delete p_cell_killer;
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
#endif /*TESTMAKENICECRYPTSIMSALEXW_HPP_*/
