#ifndef TESTCELLKILLERS_HPP_
#define TESTCELLKILLERS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "TissueSimulation.hpp"
#include "TissueCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "CellsGenerator.hpp"

class TestCellKillers : public CxxTest::TestSuite
{
public:
    void TestRandomCellKiller(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        
        // read in mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, mesh);
        
        Crypt<2> crypt(mesh, cells);
        
        // Get a reference to the cells held in crypt
        std::list<TissueCell>& r_cells = crypt.rGetCells();
        
        // bad probabilities passed in
        TS_ASSERT_THROWS_ANYTHING(RandomCellKiller<2> random_cell_killer(&crypt, -0.1));
        TS_ASSERT_THROWS_ANYTHING(RandomCellKiller<2> random_cell_killer(&crypt,  1.1));
        
        RandomCellKiller<2> random_cell_killer(&crypt, 0.05);
       
        // check that a single cell reaches apoptosis
        unsigned max_tries=0;
        while (!r_cells.begin()->HasApoptosisBegun() && max_tries<99)
        {
            random_cell_killer.TestAndLabelSingleCellForApoptosis(*r_cells.begin());
            max_tries++;
        }
        TS_ASSERT_DIFFERS(max_tries, 99u);
        TS_ASSERT_DIFFERS(max_tries, 0u);
        
        
        // check that some of the vector of cells reach apotosis
        random_cell_killer.TestAndLabelCellsForApoptosisOrDeath();
        
        std::set< double > old_locations;
        
        bool apoptosis_cell_found=false;
        std::list<TissueCell>::iterator cell_it = r_cells.begin();
        ++cell_it;
        while (cell_it != r_cells.end() && !apoptosis_cell_found)
        {
            if (cell_it->HasApoptosisBegun())
            {
                apoptosis_cell_found = true;
            }
            ++cell_it;
        }
        
        TS_ASSERT(apoptosis_cell_found);
        
        // increment time to a time after death 
        double death_time = p_simulation_time->GetDimensionalisedTime() + p_params->GetApoptosisTime();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(death_time+1.0, 1);
        p_simulation_time->IncrementTimeOneStep();
        
        // store 'locations' of cells which are not dead
        for (std::list<TissueCell>::iterator it = r_cells.begin();
             it != r_cells.end(); ++it)
        {
            if (!it->IsDead())
            {
                Node<2>* p_node = mesh.GetNode(it->GetNodeIndex());
                c_vector< double, 2 > location = p_node->rGetLocation();
                old_locations.insert(location[0]+location[1]*1000);
            }
        }
        
        // remove dead cells...
        crypt.RemoveDeadCells();
        
        // check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<TissueCell>::iterator it = r_cells.begin();
             it != r_cells.end(); ++it)
        {
            TS_ASSERT(!it->IsDead());
            Node<2>* p_node = mesh.GetNode(it->GetNodeIndex());
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
        p_params->Reset();
        
        // read in mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, mesh);
        Crypt<2> crypt(mesh, cells);
        
        p_params->SetCryptWidth(0.5);
        p_params->SetCryptLength(0.5);

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        sloughing_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

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
  
        crypt.RemoveDeadCells();

        for(Crypt<2>::Iterator iter = crypt.Begin();
            iter!=crypt.End();
            ++iter)
        {
            double x = iter.rGetLocation()[0];
            double y = iter.rGetLocation()[1];
            
            TS_ASSERT_LESS_THAN_EQUALS(x, 0.5);
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.5);
        }

        SimulationTime::Destroy();
    }   


    void TestSloughingCellKillerTopOnly(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
 
        // read in mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        std::vector<TissueCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = 0.0;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        Crypt<2> crypt(mesh, cells);
        
        p_params->SetCryptWidth(0.5);
        p_params->SetCryptLength(0.5);

        SloughingCellKiller sloughing_cell_killer(&crypt);
        sloughing_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

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
        
        crypt.RemoveDeadCells();

        for(Crypt<2>::Iterator iter = crypt.Begin();
            iter!=crypt.End();
            ++iter)
        {
            double y = iter.rGetLocation()[1];
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.5);
        }

        SimulationTime::Destroy();
    }   
    
    void TestArchivingOfRandomCellKiller() throw (Exception)
    {
        CancerParameters::Instance()->Reset();    
    
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "random_killer.arch";

        {
            // Create an ouput archive
            RandomCellKiller<2> cell_killer(NULL, 0.134);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            RandomCellKiller<2> * const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_DELTA(p_cell_killer->GetDeathProbability(), 0.134, 1e-9);
       }

       {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            RandomCellKiller<2>* p_cell_killer;

            // restore from the archive
            input_arch >> p_cell_killer;
           
            // test we have restored the probability correctly.
            TS_ASSERT_DELTA(p_cell_killer->GetDeathProbability(), 0.134, 1e-9);
            delete p_cell_killer;
        }
    }
        

    void TestArchivingOfSloughingCellKiller() throw (Exception)
    {
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "sloughing_killer.arch";

        CancerParameters *p_params = CancerParameters::Instance();
       
        p_params->SetCryptLength(10.0);
        p_params->SetCryptWidth(5.0);

        {
            // Create an ouput archive
            SloughingCellKiller cell_killer(NULL, true);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            SloughingCellKiller * const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_EQUALS(p_cell_killer->GetSloughSides(), true);
            TS_ASSERT_DELTA(p_cell_killer->GetCryptLength(), 10.0, 1e-9);
            TS_ASSERT_DELTA(p_cell_killer->GetCryptWidth(), 5.0, 1e-9);
        }
 
        // Change the cancer parameters
        p_params->SetCryptLength(12.0);
        p_params->SetCryptWidth(6.0);

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
           
            SloughingCellKiller* p_cell_killer;

            // restore from the archive
            input_arch >> p_cell_killer;
           
            // test we have restored the sloughing properties correctly.
            TS_ASSERT_EQUALS(p_cell_killer->GetSloughSides(), true);
            TS_ASSERT_DELTA(p_cell_killer->GetCryptLength(), 10.0, 1e-9);
            TS_ASSERT_DELTA(p_cell_killer->GetCryptWidth(), 5.0, 1e-9);
            delete p_cell_killer;
        }  
    }
};

#endif /*TESTCELLKILLERS_HPP_*/
