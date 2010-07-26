/*

Copyright (C) University of Oxford, 2005-2010

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

#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "CellwiseData.hpp"
#include "TrianglesMeshReader.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractCellKiller.
 */
class TestCellKillers : public AbstractCellBasedTestSuite
{
public:

    void TestRandomCellKiller() throw(Exception)
    {
        // Set up singleton classes
        TissueConfig* p_params = TissueConfig::Instance();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(32.0, 32);

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Get a reference to the cells held in tissue
        std::list<TissueCellPtr>& r_cells = tissue.rGetCells();

        // Check for bad probabilities being passed in
        TS_ASSERT_THROWS_THIS(RandomCellKiller<2> random_cell_killer(&tissue, -0.1),
                              "Probability of death must be between zero and one");

        TS_ASSERT_THROWS_THIS(RandomCellKiller<2> random_cell_killer(&tissue,  1.1),
                              "Probability of death must be between zero and one");

        // Create cell killer
        RandomCellKiller<2> random_cell_killer(&tissue, 0.05);

        // Check that a single cell reaches apoptosis
        unsigned max_tries=0;
        while (!(*r_cells.begin())->HasApoptosisBegun() && max_tries<99)
        {
            random_cell_killer.TestAndLabelSingleCellForApoptosis(*r_cells.begin());
            max_tries++;
        }
        TS_ASSERT_DIFFERS(max_tries, 99u);
        TS_ASSERT_DIFFERS(max_tries, 0u);

        // Check that some of the vector of cells reach apotosis
        random_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        std::set< double > old_locations;

        bool apoptosis_cell_found = false;
        std::list<TissueCellPtr>::iterator cell_it = r_cells.begin();
        ++cell_it;
        while (cell_it != r_cells.end() && !apoptosis_cell_found)
        {
            if ((*cell_it)->HasApoptosisBegun())
            {
                apoptosis_cell_found = true;
            }
            ++cell_it;
        }

        TS_ASSERT_EQUALS(apoptosis_cell_found, true);

        // Increment time to a time after cell death
        double death_time = p_simulation_time->GetTime() + p_params->GetApoptosisTime();
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(death_time+1.0, 1);
        p_simulation_time->IncrementTimeOneStep();

        // Store 'locations' of cells which are not dead
        for (std::list<TissueCellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            if (!(*cell_iter)->IsDead())
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
        for (std::list<TissueCellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS((*cell_iter)->IsDead(), false);
            Node<2>* p_node = tissue.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
    }


    void TestSloughingCellKillerTopAndSides() throw(Exception)
    {
        // Set up singleton classes
        TissueConfig* p_params = TissueConfig::Instance();

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);

        // Create cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        p_params->SetCryptWidth(0.5);
        p_params->SetCryptLength(0.5);

        // Create cell killer and kill cells
        SloughingCellKiller<2> sloughing_cell_killer(&tissue, true);
        sloughing_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            double x = tissue.GetLocationOfCellCentre(*cell_iter)[0];
            double y = tissue.GetLocationOfCellCentre(*cell_iter)[1];

            if ( (x<0) || (x>0.5) || (y>0.5))
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        tissue.RemoveDeadCells();

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            double x = tissue.GetLocationOfCellCentre(*cell_iter)[0];
            double y = tissue.GetLocationOfCellCentre(*cell_iter)[1];

            TS_ASSERT_LESS_THAN_EQUALS(x, 0.5);
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.5);
        }
    }


    void TestSloughingCellKillerTopOnly() throw(Exception)
    {
        // Set up singleton classes
        TissueConfig* p_params = TissueConfig::Instance();

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);

        // Create cells
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            p_cell->SetBirthTime(0.0);

            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        p_params->SetCryptWidth(0.5);
        p_params->SetCryptLength(0.5);

        // Create cell killer and kill cells
        SloughingCellKiller<2> sloughing_cell_killer(&tissue);
        sloughing_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            double y = tissue.GetLocationOfCellCentre(*cell_iter)[1];
            if (y>0.5)
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        tissue.RemoveDeadCells();

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            double y = tissue.GetLocationOfCellCentre(*cell_iter)[1];
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.5);
        }
    }


    void TestSloughingCellKillerIn1d() throw(Exception)
    {
        // Set up singleton classes
        TissueConfig* p_params = TissueConfig::Instance();

        // Create 1D mesh
        unsigned num_cells = 14;
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(num_cells-1);

        // Create cells
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            p_cell->SetBirthTime(0.0);

            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> tissue(mesh, cells);

        // Set the crypt length so that 2 cells should be sloughed off
        double crypt_length = 12.5;
        p_params->SetCryptLength(crypt_length);

        // Create cell killer and kill cells
        SloughingCellKiller<1> sloughing_cell_killer(&tissue);
        sloughing_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractTissue<1>::Iterator cell_iter = tissue.Begin();
            cell_iter != tissue.End();
            ++cell_iter)
        {
            double x = tissue.GetLocationOfCellCentre(*cell_iter)[0];
            if (x > crypt_length)
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        // Check that dead cells were correctly removed
        tissue.RemoveDeadCells();

        for (AbstractTissue<1>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            double x = tissue.GetLocationOfCellCentre(*cell_iter)[0];
            TS_ASSERT_LESS_THAN_EQUALS(x, crypt_length);
        }
    }


    void TestSloughingCellKillerIn3d() throw(Exception)
    {
        // Create 3D mesh
        MutableMesh<3,3> mesh;
        mesh.ConstructCuboid(4, 5, 6);

        // Create cells
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            p_cell->SetBirthTime(0.0);

            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<3> tissue(mesh, cells);

        // Create cell killer
        SloughingCellKiller<3> sloughing_cell_killer(&tissue);

        // Check that an exception is thrown, as this method is not yet implemented in 3D
        TS_ASSERT_THROWS_THIS(sloughing_cell_killer.TestAndLabelCellsForApoptosisOrDeath(), "SloughingCellKiller is not yet implemented in 3D");
    }


    void TestOxygenBasedCellKiller() throw(Exception)
    {
        // Set up
        TissueConfig::Instance()->SetHepaOneParameters();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 1.0;
        unsigned num_timesteps = 100*(unsigned)end_time; // ensure the time step is not too small
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Before we can do anything with the cell killer, we need to set up CellwiseData
        std::vector<double> oxygen_concentration;

        // Set the oxygen concentration to be zero
        oxygen_concentration.push_back(0.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        OxygenBasedCellKiller<2> bad_cell_killer(&tissue);

        // Get a reference to the cells held in tissue
        std::list<TissueCellPtr>& r_cells = tissue.rGetCells();

        // Reset cell types to STEM
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            cell_iter->GetCellCycleModel()->SetCellProliferativeType(STEM);
        }

        TS_ASSERT_THROWS_NOTHING(OxygenBasedCellKiller<2> oxygen_based_cell_killer(&tissue));

        OxygenBasedCellKiller<2> oxygen_based_cell_killer(&tissue);

        TS_ASSERT_THROWS_NOTHING(oxygen_based_cell_killer.TestAndLabelSingleCellForApoptosis(*r_cells.begin()));

        // Check that a single cell reaches apoptosis
        TS_ASSERT_EQUALS((*r_cells.begin())->HasApoptosisBegun(), false);

        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());
        (*r_cells.begin())->AddCellProperty(p_apoptotic_state);
        oxygen_based_cell_killer.TestAndLabelSingleCellForApoptosis(*r_cells.begin());

        TS_ASSERT((*r_cells.begin())->HasApoptosisBegun());

        // Increment time to a time after death
        p_simulation_time->IncrementTimeOneStep();

        // Store 'locations' of cells which are not dead
        std::set< double > old_locations;
        for (std::list<TissueCellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            if (!(*cell_iter)->IsDead())
            {
                Node<2>* p_node = tissue.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove the dead cell
        tissue.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<TissueCellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS((*cell_iter)->IsDead(), false);
            Node<2>* p_node = tissue.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
        CellwiseData<2>::Destroy();
    }


    void TestArchivingOfRandomCellKiller() throw (Exception)
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "random_killer.arch";

        {
            // Create an output archive
            RandomCellKiller<2> cell_killer(NULL, 0.134);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            RandomCellKiller<2>* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_DELTA(p_cell_killer->GetDeathProbabilityInAnHour(), 0.134, 1e-9);
       }

       {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            RandomCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the probability correctly
            TS_ASSERT_DELTA(p_cell_killer->GetDeathProbabilityInAnHour(), 0.134, 1e-9);

            delete p_cell_killer;
        }
    }


    void TestArchivingOfSloughingCellKiller() throw (Exception)
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "sloughing_killer.arch";

        TissueConfig* p_params = TissueConfig::Instance();

        p_params->SetCryptLength(10.0);
        p_params->SetCryptWidth(5.0);

        {
            // Create an output archive
            SloughingCellKiller<2> cell_killer(NULL, true);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            SloughingCellKiller<2>* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_EQUALS(p_cell_killer->GetSloughSides(), true);
        }

        // Change the model parameters
        p_params->SetCryptLength(12.0);
        p_params->SetCryptWidth(6.0);

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SloughingCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the sloughing properties correctly
            TS_ASSERT_EQUALS(p_cell_killer->GetSloughSides(), true);

            delete p_cell_killer;
        }
    }


    void TestArchivingOfOxygenBasedCellKiller() throw (Exception)
    {
        // Set up
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "oxygen_based_killer.arch";

        {
            // Create an output archive
            OxygenBasedCellKiller<2> cell_killer(NULL);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            OxygenBasedCellKiller<2>* const p_cell_killer = &cell_killer;

            p_cell_killer->SetHypoxicConcentration(0.3);

            output_arch << p_cell_killer;

            TS_ASSERT_DELTA(p_cell_killer->GetHypoxicConcentration(), 0.3, 1e-5);
       }

       {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            OxygenBasedCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the sloughing properties correctly
            TS_ASSERT_DELTA(p_cell_killer->GetHypoxicConcentration(), 0.3, 1e-5);

            delete p_cell_killer;
        }
    }
};

#endif /*TESTCELLKILLERS_HPP_*/
