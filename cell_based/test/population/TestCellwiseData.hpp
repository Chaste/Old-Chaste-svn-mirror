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
#ifndef TESTCELLWISEDATA_HPP_
#define TESTCELLWISEDATA_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "CellsGenerator.hpp"
#include "CellwiseData.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellLabel.hpp"

/**
 * This class contains tests for methods on the class CellwiseData.
 */
class TestCellwiseData : public AbstractCellBasedTestSuite
{
public:

    void TestCellwiseDataSimple() throw(Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);

        // One variable tests

        CellwiseData<2>* p_data = CellwiseData<2>::Instance();

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);
        TS_ASSERT_THROWS_THIS(p_data->SetCellPopulation(&cell_population),"SetCellPopulation must be called after SetNumCellsAndVars()");

        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);

        p_data->SetCellPopulation(&cell_population);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->GetNumVariables(), 1u);

        p_data->SetValue(1.23, mesh.GetNode(0)->GetIndex());
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), 1.23, 1e-12);

        p_data->SetValue(2.23, mesh.GetNode(1)->GetIndex());
        ++cell_iter;
        TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), 2.23, 1e-12);

        // Test ReallocateMemory method
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);

        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetBirthTime(-1);

        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 0.2;
        new_cell_location[1] = 0.3;
        cell_population.AddCell(p_new_cell, new_cell_location, cells[0] /*random choice of parent*/);

        TS_ASSERT_THROWS_NOTHING(p_data->ReallocateMemory());
        TS_ASSERT_EQUALS(p_data->mData.size(), cell_population.rGetMesh().GetNumNodes());

        // Coverage
        std::vector<double> constant_value;
        constant_value.push_back(1.579);
        p_data->SetConstantDataForTesting(constant_value);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), 1.579, 1e-12);
            TS_ASSERT_THROWS_THIS(p_data->GetValue(*cell_iter,1),
            		"Request for variable above mNumberOfVariables. Call SetNumCellsAndVars() to increase it.");
        }

        TS_ASSERT_THROWS_THIS(p_data->SetValue(0.0, 0, 1),
        		"Request for variable above mNumberOfVariables. Call SetNumCellsAndVars() to increase it.");

        p_data->Destroy();

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);

        // Two variable tests

        p_data = CellwiseData<2>::Instance();

        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 2);
        p_data->SetCellPopulation(&cell_population);

        TS_ASSERT_THROWS_THIS(p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1),"SetNumCellsAndVars() must be called before setting the CellPopulation (and after a Destroy)");
        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->GetNumVariables(), 2u);

        p_data->SetValue(3.23, mesh.GetNode(0)->GetIndex(), 1);
        AbstractCellPopulation<2>::Iterator cell_iter2 = cell_population.Begin();

        TS_ASSERT_DELTA(p_data->GetValue(*cell_iter2, 1), 3.23, 1e-12);

        p_data->SetValue(4.23, mesh.GetNode(1)->GetIndex(), 1);
        ++cell_iter2;

        TS_ASSERT_DELTA(p_data->GetValue(*cell_iter2, 1), 4.23, 1e-12);

        // Other values should have been initialised to zero
        ++cell_iter2;
        TS_ASSERT_DELTA(p_data->GetValue(*cell_iter2, 0), 0.0, 1e-12);

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestArchiveCellwiseData()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index, so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Work out where to put the archive
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "cellwise_data.arch";
        ArchiveLocationInfo::SetMeshFilename("cellwise_data_mesh");

        {
            // Set up the data store
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();
            p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
            p_data->SetCellPopulation(&cell_population);

            // Put some data in
            unsigned i = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                p_data->SetValue((double) i, cell_population.GetLocationIndexUsingCell(*cell_iter), 0);
                i++;
            }

            TS_ASSERT_EQUALS(p_data->IsSetUp(), true);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write to the archive
            (*p_arch) << static_cast<const CellwiseData<2>&>(*CellwiseData<2>::Instance());

            CellwiseData<2>::Destroy();
        }

        {
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> *p_data;

            // Check the data
            TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);
            TS_ASSERT_EQUALS(p_data->IsSetUp(), true);
            TS_ASSERT_EQUALS(p_data->mUseConstantDataForTesting, false);
            TS_ASSERT_EQUALS(p_data->GetNumVariables(), 1u);

            // We will have constructed a new cell population on load, so use the new cell population
            AbstractCellPopulation<2>& cell_population = p_data->rGetCellPopulation();

            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(p_data->GetValue(*cell_iter, 0u), (double) cell_population.GetLocationIndexUsingCell(*cell_iter), 1e-12);
            }

            // Tidy up
            CellwiseData<2>::Destroy();
            delete (&cell_population);
        }
    }

    void TestArchiveCellwiseDataWithVertexBasedCellPopulation()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "vertex_cellwise.arch";
        // The following line is required because the loading of a cell population
        // is usually called by the method CellBasedSimulation::Load()
        ArchiveLocationInfo::SetMeshFilename("vertex_cellwise");

        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Need to set up time
        unsigned num_steps = 10;
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2>* const p_cell_population = new VertexBasedCellPopulation<2>(mesh, cells);

        // Cells have been given birth times of 0 and -1.
        // Loop over them to run to time 0.0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
             cell_iter != p_cell_population->End();
             ++cell_iter)
        {
            cell_iter->ReadyToDivide();
        }

        {
            // Set up the data store
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();
            p_data->SetNumCellsAndVars(p_cell_population->GetNumRealCells(), 1);
            p_data->SetCellPopulation(p_cell_population);

            // Put some data in
            unsigned i = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                p_data->SetValue((double) i, p_cell_population->GetLocationIndexUsingCell(*cell_iter), 0);
                i++;
            }

            TS_ASSERT_EQUALS(p_data->IsSetUp(), true);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write to the archive
            (*p_arch) << static_cast<const CellwiseData<2>&>(*CellwiseData<2>::Instance());

            CellwiseData<2>::Destroy();
            delete p_cell_population;
        }

        {
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> *p_data;

            // Check the data
            TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);
            TS_ASSERT_EQUALS(p_data->IsSetUp(), true);

            // We will have constructed a new cell population on load, so use the new cell population
            AbstractCellPopulation<2>& cell_population = p_data->rGetCellPopulation();

            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(p_data->GetValue(*cell_iter, 0u), (double) cell_population.GetLocationIndexUsingCell(*cell_iter), 1e-12);
            }

            // Tidy up
            CellwiseData<2>::Destroy();
            delete (&cell_population);
        }
    }
};

#endif /*TESTCELLWISEDATA_HPP_*/
