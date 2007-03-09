#ifndef TESTCELLKILLERS_HPP_
#define TESTCELLKILLERS_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "CryptSimulation2DPeriodic.hpp"
#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "CryptHoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"
#include "RandomCellKiller.hpp"


class TestCellKillers : public CxxTest::TestSuite
{
public:
    void TestRandomCellKiller(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen; // passed into crypt sim for coverage
        
        // read in mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells by iterating through the mesh nodes
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // Any old rubbish here just so the simulation time is set up for cell cycle models
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 9);
        
        
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            
            double y = mesh.GetNode(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -random_num_gen.ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 6.5)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 8)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = -1; //hours
            }
            
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(num_cells-i-1);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
            
            
        }
        
        RandomCellKiller<2> random_cell_killer;
        random_cell_killer.SetCellsAndMesh(&cells, &mesh);
        
        // check that a single cell reaches apoptosis    
        unsigned max_tries=0;
        while(!cells[0].HasApoptosisBegun() && max_tries<99)
        {
            random_cell_killer.TestAndLabelSingleCellForApoptosis(0);
            max_tries++;
        }
        TS_ASSERT_DIFFERS(max_tries, 99u);
        TS_ASSERT_DIFFERS(max_tries, 0u);
        
        
        // check that some of the vector of cells reach apotosis
        random_cell_killer.TestAndLabelCellsForApoptosis();
        
        std::set< double > old_locations;
        
        bool apoptosis_cell_found=false;
        unsigned count=1;
        while(count<cells.size() && apoptosis_cell_found ==  false)
        {
            if (cells[count].HasApoptosisBegun())
            {
                apoptosis_cell_found = true;
            }
            count++;
        }
            
        TS_ASSERT(apoptosis_cell_found);


        double death_time = p_simulation_time->GetDimensionalisedTime() + p_params->GetApoptosisTime();
        while (p_simulation_time->GetDimensionalisedTime() < death_time)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        
        // store 'locations' of cells which are not dead
        for (count=0; count<cells.size(); count++)
        {
            if (!cells[count].IsDead())
            {
                Node<2>* p_node = mesh.GetNode(cells[count].GetNodeIndex());
                c_vector< double, 2 > location = p_node->rGetLocation();
                old_locations.insert(location[0]+location[1]*1000);
            }  
        }
        
        
        // remove dead cells...
        random_cell_killer.RemoveDeadCells();
        
        // check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for( count=0; count<cells.size(); count++)
        {
            TS_ASSERT(!cells[count].IsDead());
            Node<2>* p_node = mesh.GetNode(cells[count].GetNodeIndex());
            c_vector< double, 2 > location = p_node->rGetLocation();
            new_locations.insert(location[0]+location[1]*1000);
        }
        
        TS_ASSERT(new_locations == old_locations);
        
   }
    
    
};

#endif /*TESTCELLKILLERS_HPP_*/
