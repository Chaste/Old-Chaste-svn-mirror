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
#ifndef TESTCELLWISEDATA_HPP_
#define TESTCELLWISEDATA_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "CellsGenerator.hpp"
#include "CellwiseData.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "WildTypeCellMutationState.hpp"

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
        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);

        // One variable tests

        CellwiseData<2>* p_data = CellwiseData<2>::Instance();

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);
        TS_ASSERT_THROWS_THIS(p_data->SetTissue(&tissue),"SetTissue must be called after SetNumCellsAndVars()");

        p_data->SetNumCellsAndVars(tissue.GetNumRealCells(), 1);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);

        p_data->SetTissue(&tissue);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->GetNumVariables(), 1u);

        p_data->SetValue(1.23, mesh.GetNode(0)->GetIndex());
        AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
        TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), 1.23, 1e-12);

        p_data->SetValue(2.23, mesh.GetNode(1)->GetIndex());
        ++cell_iter;
        TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), 2.23, 1e-12);

        // Test ReallocateMemory method
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);

        TissueCell new_cell(p_state, p_model);
        new_cell.SetBirthTime(-1);

        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 0.2;
        new_cell_location[1] = 0.3;
        tissue.AddCell(new_cell, new_cell_location);

        TS_ASSERT_THROWS_NOTHING(p_data->ReallocateMemory());
        TS_ASSERT_EQUALS(p_data->mData.size(), tissue.rGetMesh().GetNumNodes());

        // Coverage
        std::vector<double> constant_value;
        constant_value.push_back(1.579);
        p_data->SetConstantDataForTesting(constant_value);

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
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

        p_data->SetNumCellsAndVars(tissue.GetNumRealCells(), 2);
        p_data->SetTissue(&tissue);

        TS_ASSERT_THROWS_THIS(p_data->SetNumCellsAndVars(tissue.GetNumRealCells(), 1),"SetNumCellsAndVars() must be called before setting the Tissue (and after a Destroy)");
        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->GetNumVariables(), 2u);

        p_data->SetValue(3.23, mesh.GetNode(0)->GetIndex(), 1);
        AbstractTissue<2>::Iterator cell_iter2 = tissue.Begin();

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
        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a tissue
        MeshBasedTissue<2> tissue(mesh,cells);

        // Work out where to put the archive
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "cellwise_data.arch";
        ArchiveLocationInfo::SetMeshFilename("cellwise_data_mesh");

        {
            // Set up the data store
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();
            p_data->SetNumCellsAndVars(tissue.GetNumRealCells(), 1);
            p_data->SetTissue(&tissue);

            // Put some data in
            unsigned i=0;
            for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
                 cell_iter != tissue.End();
                 ++cell_iter)
            {
                p_data->SetValue((double) i, tissue.GetLocationIndexUsingCell(*cell_iter), 0);
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

            // We will have constructed a new tissue on load, so use the new tissue
            AbstractTissue<2>& tissue = p_data->rGetTissue();

            for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
                 cell_iter != tissue.End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(p_data->GetValue(*cell_iter, 0u), (double) tissue.GetLocationIndexUsingCell(*cell_iter), 1e-12);
            }

            // Tidy up
            CellwiseData<2>::Destroy();
            delete (&tissue);
        }
    }

};

#endif /*TESTCELLWISEDATA_HPP_*/
