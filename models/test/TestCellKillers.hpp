#ifndef TESTCELLKILLERS_HPP_
#define TESTCELLKILLERS_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "TissueSimulation.hpp"
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
#include "SloughingCellKiller.hpp"

class TestCellKillers : public CxxTest::TestSuite
{
public:
    void TestRandomCellKiller(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        
        // read in mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = 0.0;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        Crypt<2> crypt(mesh, cells);
        // Get a reference to the cells held in crypt
        std::vector<MeinekeCryptCell>& r_cells = crypt.rGetCells();
        
        // bad probabilities passed in
        TS_ASSERT_THROWS_ANYTHING(RandomCellKiller<2> random_cell_killer(&crypt, -0.1));
        TS_ASSERT_THROWS_ANYTHING(RandomCellKiller<2> random_cell_killer(&crypt,  1.1));
        
        RandomCellKiller<2> random_cell_killer(&crypt, 0.05);
       
        // check that a single cell reaches apoptosis
        unsigned max_tries=0;
        while (!r_cells[0].HasApoptosisBegun() && max_tries<99)
        {
            random_cell_killer.TestAndLabelSingleCellForApoptosis(r_cells[0]);
            max_tries++;
        }
        TS_ASSERT_DIFFERS(max_tries, 99u);
        TS_ASSERT_DIFFERS(max_tries, 0u);
        
        
        // check that some of the vector of cells reach apotosis
        random_cell_killer.TestAndLabelCellsForApoptosis();
        
        std::set< double > old_locations;
        
        bool apoptosis_cell_found=false;
        unsigned count=1;
        while (count<cells.size() && apoptosis_cell_found ==  false)
        {
            if (r_cells[count].HasApoptosisBegun())
            {
                apoptosis_cell_found = true;
            }
            count++;
        }
        
        TS_ASSERT(apoptosis_cell_found);
        
        // increment time to a time after death 
        double death_time = p_simulation_time->GetDimensionalisedTime() + p_params->GetApoptosisTime();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(death_time+1.0, 1);
        p_simulation_time->IncrementTimeOneStep();
        
        // store 'locations' of cells which are not dead
        for (count=0; count<r_cells.size(); count++)
        {
            if (!r_cells[count].IsDead())
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
        for ( count=0; count<r_cells.size(); count++)
        {
            TS_ASSERT(!r_cells[count].IsDead());
            Node<2>* p_node = mesh.GetNode(r_cells[count].GetNodeIndex());
            c_vector< double, 2 > location = p_node->rGetLocation();
            new_locations.insert(location[0]+location[1]*1000);
        }
        
        TS_ASSERT(new_locations == old_locations);
        RandomNumberGenerator::Destroy();  
        SimulationTime::Destroy();
    }   

    void TestSloughingCellKillerTopAndSides(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        
        // read in mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);

        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = 0.0;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        Crypt<2> crypt(mesh, cells);
        
        p_params->SetCryptWidth(0.5);
        p_params->SetCryptLength(0.5);

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        sloughing_cell_killer.TestAndLabelCellsForApoptosis();

        for(Crypt<2>::Iterator iter = crypt.Begin();
            iter!=crypt.End();
            ++iter)
        {
            double x = iter.rGetLocation()[0];
            double y = iter.rGetLocation()[1];
            
            if( (x<0) || (x>0.5) || (y>0.5))
            {
                TS_ASSERT_EQUALS(iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(iter->IsDead(), false);
            }
        }
  
///\todo: fix this. get 'failure to delete node' error, prob because 
// trying to delete too many at once?
        
//        sloughing_cell_killer.RemoveDeadCells();
//
//        for(Crypt<2>::Iterator iter = crypt.Begin();
//            iter!=crypt.End();
//            ++iter)
//        {
//            double x = iter.rGetLocation()[0];
//            double y = iter.rGetLocation()[1];
//            
//            TS_ASSERT_LESS_THAN(x, 0.5);
//            TS_ASSERT_LESS_THAN(y, 0.5);
//        }

        SimulationTime::Destroy();
    }   


    void TestSloughingCellKillerTopOnly(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        
        // read in mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = 0.0;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        Crypt<2> crypt(mesh, cells);
        
        p_params->SetCryptWidth(0.5);
        p_params->SetCryptLength(0.5);

        SloughingCellKiller sloughing_cell_killer(&crypt);
        sloughing_cell_killer.TestAndLabelCellsForApoptosis();

        for(Crypt<2>::Iterator iter = crypt.Begin();
            iter!=crypt.End();
            ++iter)
        {
            double y = iter.rGetLocation()[1];
            if(y>0.5)
            {
                TS_ASSERT_EQUALS(iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(iter->IsDead(), false);
            }
        }
  
///\todo: fix this. get 'failure to delete node' error, prob because 
// trying to delete too many at once?
        
//        sloughing_cell_killer.RemoveDeadCells();
//
//        for(Crypt<2>::Iterator iter = crypt.Begin();
//            iter!=crypt.End();
//            ++iter)
//        {
//            double y = iter.rGetLocation()[1];
//            TS_ASSERT_LESS_THAN(y, 0.5);
//        }

        SimulationTime::Destroy();
    }   

};

#endif /*TESTCELLKILLERS_HPP_*/
