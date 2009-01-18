/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef TESTCELLKILLERS_HPP_
#define TESTCELLKILLERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "FixedCellCycleModelCellsGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "CellwiseData.hpp"
#include "TrianglesMeshReader.hpp"
#include "AbstractCancerTestSuite.hpp"

/**
 * This class contains tests for methods on classes 
 * inheriting from AbstractCellKiller.
 */
class TestCellKillers : public AbstractCancerTestSuite
{
public:

    void TestRandomCellKiller(void) throw(Exception)
    {
        // Set up
        CancerParameters *p_params = CancerParameters::Instance();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh);

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Get a reference to the cells held in tissue
        std::list<TissueCell>& r_cells = tissue.rGetCells();

        // Bad probabilities passed in
        TS_ASSERT_THROWS_ANYTHING(RandomCellKiller<2> random_cell_killer(&tissue, -0.1));
        TS_ASSERT_THROWS_ANYTHING(RandomCellKiller<2> random_cell_killer(&tissue,  1.1));

        // Create cell killer 
        RandomCellKiller<2> random_cell_killer(&tissue, 0.05);

        // Check that a single cell reaches apoptosis
        unsigned max_tries=0;
        while (!r_cells.begin()->HasApoptosisBegun() && max_tries<99)
        {
            random_cell_killer.TestAndLabelSingleCellForApoptosis(*r_cells.begin());
            max_tries++;
        }
        TS_ASSERT_DIFFERS(max_tries, 99u);
        TS_ASSERT_DIFFERS(max_tries, 0u);

        // Check that some of the vector of cells reach apotosis
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

        // Increment time to a time after cell death
        double death_time = p_simulation_time->GetTime() + p_params->GetApoptosisTime();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(death_time+1.0, 1);
        p_simulation_time->IncrementTimeOneStep();

        // Store 'locations' of cells which are not dead
        for (std::list<TissueCell>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            if (!cell_iter->IsDead())
            {
                Node<2>* p_node = tissue.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove dead cells
        tissue.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<TissueCell>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            TS_ASSERT(!cell_iter->IsDead());
            Node<2>* p_node = tissue.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
    }


    void TestSloughingCellKillerTopAndSides(void) throw(Exception)
    {
        // Set up
        CancerParameters *p_params = CancerParameters::Instance();

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);
        
        // Create cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh);
        
        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        p_params->SetCryptWidth(0.5);
        p_params->SetCryptLength(0.5);

        // Create cell killer and kill cells
        SloughingCellKiller sloughing_cell_killer(&tissue, true);
        sloughing_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractTissue<2>::Iterator iter=tissue.Begin();
            iter!=tissue.End();
            ++iter)
        {
            double x = tissue.GetNodeCorrespondingToCell(*iter)->rGetLocation()[0];
            double y = tissue.GetNodeCorrespondingToCell(*iter)->rGetLocation()[1];

            if ( (x<0) || (x>0.5) || (y>0.5))
            {
                TS_ASSERT_EQUALS(iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(iter->IsDead(), false);
            }
        }

        tissue.RemoveDeadCells();

        for (AbstractTissue<2>::Iterator iter=tissue.Begin();
            iter!=tissue.End();
            ++iter)
        {
            double x = tissue.GetNodeCorrespondingToCell(*iter)->rGetLocation()[0];
            double y = tissue.GetNodeCorrespondingToCell(*iter)->rGetLocation()[1];

            TS_ASSERT_LESS_THAN_EQUALS(x, 0.5);
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.5);
        }
    }


    void TestSloughingCellKillerTopOnly(void) throw(Exception)
    {
        // Set up
        CancerParameters *p_params = CancerParameters::Instance();

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);

        // Create cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            double birth_time = 0.0;
            cell.SetLocationIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        p_params->SetCryptWidth(0.5);
        p_params->SetCryptLength(0.5);

        // Create cell killer and kill cells
        SloughingCellKiller sloughing_cell_killer(&tissue);
        sloughing_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractTissue<2>::Iterator iter = tissue.Begin();
            iter!=tissue.End();
            ++iter)
        {
            double y = tissue.GetNodeCorrespondingToCell(*iter)->rGetLocation()[1];
            if (y>0.5)
            {
                TS_ASSERT_EQUALS(iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(iter->IsDead(), false);
            }
        }

        tissue.RemoveDeadCells();

        for (AbstractTissue<2>::Iterator iter=tissue.Begin();
            iter!=tissue.End();
            ++iter)
        {
            double y = tissue.GetNodeCorrespondingToCell(*iter)->rGetLocation()[1];
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.5);
        }
    }

    void TestArchivingOfRandomCellKiller() throw (Exception)
    {
        // Set up
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "random_killer.arch";

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

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the probability correctly
            TS_ASSERT_DELTA(p_cell_killer->GetDeathProbability(), 0.134, 1e-9);
            
            delete p_cell_killer;
        }
    }


    void TestArchivingOfSloughingCellKiller() throw (Exception)
    {
        // Set up
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "sloughing_killer.arch";

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
        }

        // Change the cancer parameters
        p_params->SetCryptLength(12.0);
        p_params->SetCryptWidth(6.0);

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SloughingCellKiller* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the sloughing properties correctly
            TS_ASSERT_EQUALS(p_cell_killer->GetSloughSides(), true);
            
            delete p_cell_killer;
        }
    }

};

#endif /*TESTCELLKILLERS_HPP_*/
