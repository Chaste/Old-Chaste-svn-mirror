#ifndef TESTMAKENICECRYPTSIMS_HPP_
#define TESTMAKENICECRYPTSIMS_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "MeinekeCryptCell.hpp"
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

// Possible types of Cell Cycle Model (just for CreateVectorOfCells method)
typedef enum CellCycleType_
{
    FIXED,
    STOCHASTIC,
    WNT,
    STOCHASTIC_WNT,
    TYSONNOVAK
} CellCycleType;


class TestNiceCryptSims : public CxxTest::TestSuite
{
    
    void CreateVectorOfCells(std::vector<MeinekeCryptCell>& rCells, 
                             ConformingTetrahedralMesh<2,2>& rMesh, 
                             CellCycleType cycleType, 
                             bool randomBirthTimes,
                             double y0 = 0.3,
                             double y1 = 2.0,
                             double y2 = 3.0,
                             double y3 = 4.0)
    {
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        unsigned num_cells = rMesh.GetNumNodes();

        AbstractCellCycleModel* p_cell_cycle_model = NULL;
        double typical_transit_cycle_time;
        double typical_stem_cycle_time;
        
        CancerParameters* p_params = CancerParameters::Instance();
        
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;

            double y = rMesh.GetNode(i)->GetPoint().rGetLocation()[1];
            
            if (cycleType==FIXED)
            {
                p_cell_cycle_model = new FixedCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==STOCHASTIC)
            {
                p_cell_cycle_model = new StochasticCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==WNT)
            {
                WntGradient wnt_gradient(OFFSET_LINEAR);
                double wnt = wnt_gradient.GetWntLevel(y);
                p_cell_cycle_model = new WntCellCycleModel(wnt,0);
                typical_transit_cycle_time = 16.0;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else if (cycleType==STOCHASTIC_WNT)
            {
                WntGradient wnt_gradient(OFFSET_LINEAR);
                double wnt = wnt_gradient.GetWntLevel(y);
                p_cell_cycle_model = new StochasticWntCellCycleModel(wnt,0);
                typical_transit_cycle_time = 16.0;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else if (cycleType==TYSONNOVAK)
            {
                p_cell_cycle_model = new TysonNovakCellCycleModel();
                typical_transit_cycle_time = 1.25;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else
            {
                EXCEPTION("Cell Cycle Type is not recognised");   
            }
            
            
            double birth_time = 0.0;
            
            if (y <= y0)
            {
                cell_type = STEM;
                generation = 0;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_stem_cycle_time; // hours
                }
            }
            else if (y < y1)
            {
                cell_type = TRANSIT;
                generation = 1;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y2)
            {
                cell_type = TRANSIT;
                generation = 2;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y3)
            {
                cell_type = TRANSIT;
                generation = 3;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else
            {
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
                if(cycleType==WNT || cycleType==TYSONNOVAK)
                {
                    // There are no fully differentiated cells!
                    cell_type = TRANSIT;
                    
                }
                else
                {
                    cell_type = DIFFERENTIATED;
                }                
                generation = 4;
            }

             MeinekeCryptCell cell(cell_type, HEALTHY, generation, p_cell_cycle_model);
            
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            rCells.push_back(cell);
        }
    }
public:

void TestNiceCryptSimulationWithWntDependentBirthAndSloughingDeath() throw (Exception)
    {        
        CancerParameters* p_params = CancerParameters::Instance();
        std::string output_directory = "NiceCryptSim";
        
        unsigned cells_across = 12;
        unsigned cells_up = 25;
        double crypt_width = 12.0;
        unsigned thickness_of_ghost_layer = 3;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, STOCHASTIC_WNT, true);
              
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        TissueSimulation<2> simulator(crypt);
        simulator.SetOutputDirectory(output_directory);
                
        // Set length of simulation here
        simulator.SetEndTime(50.0);
        
        simulator.SetMaxCells(1000);
        simulator.SetMaxElements(2000);
        simulator.SetWntGradient(OFFSET_LINEAR);

        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&simulator.rGetCrypt(),0.01);
        simulator.AddCellKiller(p_cell_killer);
        
        // UNUSUAL SET UP HERE
        
        p_params->SetDampingConstantNormal(1.0);    // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());
        
        p_params->SetSpringStiffness(30.0); //normally 15.0;
        // 0.3/30 = 0.01 (i.e. Meineke's values)
        
        simulator.UseNonFlatBottomSurface();
        
        // END OF UNUSUAL SET UP!
        
        simulator.Solve();
        
        // set a cell to be labelled (a stemish cell)
        simulator.rGetCrypt().rGetCellAtNodeIndex(37).SetMutationState(LABELLED);
        simulator.rGetCrypt().rGetCellAtNodeIndex(194).SetMutationState(LABELLED);
        simulator.rGetCrypt().rGetCellAtNodeIndex(304).SetMutationState(LABELLED);
        
        simulator.Save(); 

        for (double t=50; t<450.5; t += 50)
        {
            TissueSimulation<2>* p_simulator = TissueSimulation<2>::Load("NiceCryptSim",t);
            p_simulator->SetEndTime(t+50);
            p_simulator->Solve();
            p_simulator->Save();
            delete p_simulator;
        }
                
        delete p_cell_killer;
        delete p_params;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

};

#endif /*TESTMAKENICECRYPTSIMS_HPP_*/
