/*

Copyright (C) University of Oxford, 2005-2011

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
#include "TargetedCellKiller.hpp"
#include "RandomCellKiller.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"
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

    void TestTargetedCellKiller() throw(Exception)
    {
        // Set up singleton classes
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(32.0, 32);

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        // Create cell killer
        TargetedCellKiller<2> single_cell_killer(&cell_population, 1u);

        TS_ASSERT_EQUALS(single_cell_killer.GetIdentifier(), "TargetedCellKiller-2");

        // Check that some of the vector of cells reach apotosis
        single_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        std::set<double> old_locations;

        std::list<CellPtr>::iterator cell_it = r_cells.begin();
        TS_ASSERT(!(*cell_it)->IsDead());
        ++cell_it;
        TS_ASSERT((*cell_it)->IsDead());
        ++cell_it;

        while (cell_it != r_cells.end())
        {
            TS_ASSERT(!(*cell_it)->IsDead());
            ++cell_it;
        }

        // Store 'locations' of cells which are not dead
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            if (!(*cell_iter)->IsDead())
            {
                Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove dead cells
        cell_population.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS((*cell_iter)->IsDead(), false);
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
    }

    void TestRandomCellKiller() throw(Exception)
    {
        // Set up singleton classes
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(32.0, 32);

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        double death_time = p_simulation_time->GetTime() + cells[0]->GetApoptosisTime();

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        // Check for bad probabilities being passed in
        TS_ASSERT_THROWS_THIS(RandomCellKiller<2> random_cell_killer(&cell_population, -0.1),
                             "Probability of death must be between zero and one");

        TS_ASSERT_THROWS_THIS(RandomCellKiller<2> random_cell_killer(&cell_population,  1.1),
                             "Probability of death must be between zero and one");

        // Create cell killer
        RandomCellKiller<2> random_cell_killer(&cell_population, 0.05);

        TS_ASSERT_EQUALS(random_cell_killer.GetIdentifier(), "RandomCellKiller-2");

        // Check that a single cell reaches apoptosis
        unsigned max_tries = 0;
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
        std::list<CellPtr>::iterator cell_it = r_cells.begin();
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
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(death_time+1.0, 1);
        p_simulation_time->IncrementTimeOneStep();

        // Store 'locations' of cells which are not dead
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
            cell_iter != r_cells.end();
            ++cell_iter)
        {
            if (!(*cell_iter)->IsDead())
            {
                Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove dead cells
        cell_population.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
            cell_iter != r_cells.end();
            ++cell_iter)
        {
            TS_ASSERT_EQUALS((*cell_iter)->IsDead(), false);
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
    }

    void TestOxygenBasedCellKiller() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 1.0;
        unsigned num_timesteps = 100*(unsigned)end_time; // ensure the time step is not too small
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Set some model parameters for the cell-cycle model
        for (unsigned index=0; index < cells.size(); index++)
        {
            cells[index]->GetCellCycleModel()->SetStemCellG1Duration(8.0);
            cells[index]->GetCellCycleModel()->SetTransitCellG1Duration(8.0);
        }

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Before we can do anything with the cell killer, we need to set up CellwiseData
        std::vector<double> oxygen_concentration;

        // Set the oxygen concentration to be zero
        oxygen_concentration.push_back(0.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        OxygenBasedCellKiller<2> bad_cell_killer(&cell_population);

        // Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        // Reset cell types to STEM
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellCycleModel()->SetCellProliferativeType(STEM);
        }

        TS_ASSERT_THROWS_NOTHING(OxygenBasedCellKiller<2> oxygen_based_cell_killer(&cell_population));

        OxygenBasedCellKiller<2> oxygen_based_cell_killer(&cell_population);

        TS_ASSERT_EQUALS(oxygen_based_cell_killer.GetIdentifier(), "OxygenBasedCellKiller-2");

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
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            if (!(*cell_iter)->IsDead())
            {
                Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove the dead cell
        cell_population.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS((*cell_iter)->IsDead(), false);
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
        CellwiseData<2>::Destroy();
    }

    void TestPlaneBasedCellKillerIn1d() throw(Exception)
    {
        // Create 1D mesh
        unsigned num_cells = 14;
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(num_cells-1);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<1> cell_population(mesh, cells);

        // Create cell killer and kill cells
        c_vector<double, 1> point;
        point(0) = 10.0;
        PlaneBasedCellKiller<1> cell_killer(&cell_population, point, unit_vector<double>(1,0)); // x<10
        cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<1>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            if (x > point(0))
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        // Check that dead cells were correctly removed
        cell_population.RemoveDeadCells();

        for (AbstractCellPopulation<1>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            TS_ASSERT_LESS_THAN_EQUALS(x, point(0));
        }
    }


    void TestPlaneBasedCellKillerIn2d() throw(Exception)
    {
        // Create mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(-0.25,-0.25);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Create cell killer and kill cells
        PlaneBasedCellKiller<2> cell_killer(&cell_population, zero_vector<double>(2), unit_vector<double>(2,1)); // y<0
        cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            if (y > 0.0)
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        cell_population.RemoveDeadCells();

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            TS_ASSERT_LESS_THAN_EQUALS(y, 0.0);
        }
    }

    void TestPlaneBasedCellKillerIn3d() throw(Exception)
    {
        // Create 3D mesh
        MutableMesh<3,3> mesh;
        mesh.ConstructCuboid(4, 5, 6);
        mesh.Translate(-2.0,-2.0, -2.0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        // Create cell killer
        PlaneBasedCellKiller<3> cell_killer(&cell_population,  zero_vector<double>(3), unit_vector<double>(3,2)); // z<0
        cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        // Check that cells were labelled for death correctly
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double z = cell_population.GetLocationOfCellCentre(*cell_iter)[2];
            if (z > 0.0)
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        cell_population.RemoveDeadCells();

        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double z = cell_population.GetLocationOfCellCentre(*cell_iter)[2];
            TS_ASSERT_LESS_THAN_EQUALS(z, 0.0);
        }
    }



    void TestArchivingOfTargetedCellKiller() throw (Exception)
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "single_cell_killer.arch";

        {
            // Create an output archive
             TargetedCellKiller<2> cell_killer(NULL, 1u);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            TargetedCellKiller<2>* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_EQUALS(p_cell_killer->GetTargetIndex(), 1u);
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            TargetedCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the Target Cell correctly
            TS_ASSERT_EQUALS(p_cell_killer->GetTargetIndex(), 1u);

            delete p_cell_killer;
       }
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

            output_arch << p_cell_killer;
       }

       {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            OxygenBasedCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            TS_ASSERT(p_cell_killer != NULL);

            // Tidy up
            delete p_cell_killer;
        }
    }

    void TestArchivingOfPlaneBasedCellKiller() throw (Exception)
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "region_based_killer.arch";

        {
            // Create an output archive

            PlaneBasedCellKiller<2> cell_killer(NULL, zero_vector<double>(2), unit_vector<double>(2,1));

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            PlaneBasedCellKiller<2>* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            TS_ASSERT_EQUALS(p_cell_killer->rGetPointOnPlane()[0], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetPointOnPlane()[1], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetNormalToPlane()[0], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetNormalToPlane()[1], 1.0);
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            PlaneBasedCellKiller<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the region properties correctly
            TS_ASSERT_EQUALS(p_cell_killer->rGetPointOnPlane()[0], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetPointOnPlane()[1], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetNormalToPlane()[0], 0.0);
            TS_ASSERT_EQUALS(p_cell_killer->rGetNormalToPlane()[1], 1.0);

            delete p_cell_killer;
        }
    }

    void TestCellKillersOutputParameters()
    {
        std::string output_directory = "TestCellKillersOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with TargetedCellKiller
        TargetedCellKiller<2> targeted_cell_killer(NULL, 1u);
        TS_ASSERT_EQUALS(targeted_cell_killer.GetIdentifier(), "TargetedCellKiller-2");

        out_stream targeted_cell_killer_parameter_file = output_file_handler.OpenOutputFile("targeted_results.parameters");
        targeted_cell_killer.OutputCellKillerParameters(targeted_cell_killer_parameter_file);
        targeted_cell_killer_parameter_file->close();

        std::string targeted_cell_killer_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + targeted_cell_killer_results_dir + "targeted_results.parameters cell_based/test/data/TestCellKillers/targeted_results.parameters").c_str()), 0);

        // Test with RandomCellKiller
        RandomCellKiller<2> random_cell_killer(NULL, 0.01);
        TS_ASSERT_EQUALS(random_cell_killer.GetIdentifier(), "RandomCellKiller-2");

        out_stream random_cell_killer_parameter_file = output_file_handler.OpenOutputFile("random_results.parameters");
        random_cell_killer.OutputCellKillerParameters(random_cell_killer_parameter_file);
        random_cell_killer_parameter_file->close();

        std::string random_cell_killer_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + random_cell_killer_results_dir + "random_results.parameters cell_based/test/data/TestCellKillers/random_results.parameters").c_str()), 0);

        // Test with OxygenBasedCellKiller
        OxygenBasedCellKiller<2> oxygen_cell_killer(NULL);
        TS_ASSERT_EQUALS(oxygen_cell_killer.GetIdentifier(), "OxygenBasedCellKiller-2");

        out_stream oxygen_cell_killer_parameter_file = output_file_handler.OpenOutputFile("oxygen_results.parameters");
        oxygen_cell_killer.OutputCellKillerParameters(oxygen_cell_killer_parameter_file);
        oxygen_cell_killer_parameter_file->close();

        std::string oxygen_cell_killer_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + oxygen_cell_killer_results_dir + "oxygen_results.parameters cell_based/test/data/TestCellKillers/oxygen_results.parameters").c_str()), 0);

        // Test with PlaneBasedCellKiller
        PlaneBasedCellKiller<2> region_cell_killer(NULL, zero_vector<double>(2), unit_vector<double>(2,1)); // y<0;
        TS_ASSERT_EQUALS(region_cell_killer.GetIdentifier(), "PlaneBasedCellKiller-2");

        out_stream region_cell_killer_parameter_file = output_file_handler.OpenOutputFile("region_results.parameters");
        region_cell_killer.OutputCellKillerParameters(region_cell_killer_parameter_file);
        region_cell_killer_parameter_file->close();

        std::string region_cell_killer_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + region_cell_killer_results_dir + "region_results.parameters cell_based/test/data/TestCellKillers/region_results.parameters").c_str()), 0);
    }

};

#endif /*TESTCELLKILLERS_HPP_*/