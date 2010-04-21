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
#ifndef TESTMESHBASEDTISSUE_HPP_
#define TESTMESHBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedTissue.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "LabelledCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

class TestMeshBasedTissue : public AbstractCellBasedTestSuite
{
private:

    template<unsigned DIM>
    void TestSmallMeshBasedTissue(std::string meshFilename)
    {
        // Create a simple mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader(meshFilename);
        MutableMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, DIM> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create the tissue
        unsigned num_cells = cells.size();
        MeshBasedTissue<DIM> tissue(mesh, cells);

        TS_ASSERT_EQUALS(tissue.rGetMesh().GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), num_cells);

        unsigned counter = 0;
        for (typename AbstractTissue<DIM>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(*cell_iter), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            counter++;
        }

        TS_ASSERT_EQUALS(counter, tissue.GetNumRealCells());
    }


public:

    // Test construction, accessors and Iterator
    void TestSmallMeshBasedTissue1d2d3d() throw(Exception)
    {
        TestSmallMeshBasedTissue<1>("mesh/test/data/1D_0_to_1_10_elements");
        TestSmallMeshBasedTissue<2>("mesh/test/data/square_4_elements");
        TestSmallMeshBasedTissue<3>("mesh/test/data/cube_136_elements");
    }

    void TestValidateMeshBasedTissue()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node apart from one.
        // Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<mesh.GetNumNodes()-1; i++)
        {
            AbstractCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
            TissueCell cell(STEM, p_state, p_cell_cycle_model);
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Fails as no cell corresponding to node 4
        std::vector<TissueCell> cells_copy(cells);
        TS_ASSERT_THROWS_THIS(MeshBasedTissue<2> tissue2(mesh, cells_copy),
                              "Node 4 does not appear to have a cell associated with it");

        // Add another cell
        AbstractCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
        TissueCell cell(STEM, p_state, p_cell_cycle_model);
        double birth_time = -4.0;
        cell.SetBirthTime(birth_time);
        cells.push_back(cell);

        std::vector<TissueCell> cells_copy2(cells);
        TS_ASSERT_THROWS_NOTHING(MeshBasedTissue<2> tissue2(mesh, cells_copy2));

        // A bit of Northern compatibility testing hidden here (not relevant to this test!)
        TS_ASSERT_THROWS_NOWT(MeshBasedTissue<2> tissue2(mesh, cells));
        TS_ASSERT_CHAMPION(true);
    }

    void TestCreateCellPair()
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

        // Give cells 0 and 1 specific mutations to enable later testing
        boost::shared_ptr<AbstractCellMutationState> p_labelled(new LabelledCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        cells[0].SetMutationState(p_labelled);
        cells[1].SetMutationState(p_apc1);

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Create cell pair
        std::set<TissueCell*> cell_pair = tissue.CreateCellPair(cells[0], cells[1]);

        // Check the cell pair was created correctly
        std::set<TissueCell*>::iterator cell_pair_iter = cell_pair.begin();

        TissueCell* p_cell0 = *cell_pair_iter;
        TS_ASSERT_EQUALS(p_cell0->GetMutationState(), p_labelled);

        ++cell_pair_iter;
        TissueCell* p_cell1 = *cell_pair_iter;
        TS_ASSERT_EQUALS(p_cell1->GetMutationState()->IsType<ApcOneHitCellMutationState>(), true);
    }

    void TestGetDampingConstant()
    {
        // Change the mutant damping constant to be different from the normal
        TissueConfig::Instance()->SetDampingConstantMutant(23.57);

        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Bestow mutations on some cells
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_labelled(new LabelledCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);

        cells[0].SetMutationState(p_state);
        cells[1].SetMutationState(p_apc1);
        cells[2].SetMutationState(p_apc2);
        cells[3].SetMutationState(p_bcat1);
        cells[4].SetMutationState(p_labelled);

        // Create tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        TS_ASSERT_EQUALS(tissue.UseAreaBasedDampingConstant(), false);

        double damping_const_0 = tissue.GetDampingConstant(0);
        double damping_const_1 = tissue.GetDampingConstant(1);
        double damping_const_2 = tissue.GetDampingConstant(2);
        double damping_const_3 = tissue.GetDampingConstant(3);
        double damping_const_4 = tissue.GetDampingConstant(4);

        // Check that each mutation state gives the correct damping constant
        TS_ASSERT_DELTA(damping_const_0, TissueConfig::Instance()->GetDampingConstantNormal(), 1e-6);
        TS_ASSERT_DELTA(damping_const_1, TissueConfig::Instance()->GetDampingConstantMutant(), 1e-6);
        TS_ASSERT_DELTA(damping_const_2, TissueConfig::Instance()->GetDampingConstantMutant(), 1e-6);
        TS_ASSERT_DELTA(damping_const_3, TissueConfig::Instance()->GetDampingConstantMutant(), 1e-6);
        TS_ASSERT_DELTA(damping_const_4, TissueConfig::Instance()->GetDampingConstantMutant(), 1e-6);
    }


    void TestAreaBasedDampingConstant()
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        cells[9].SetMutationState(p_apc2);

        MeshBasedTissue<2> tissue(*p_mesh, cells);

        TS_ASSERT_EQUALS(tissue.UseAreaBasedDampingConstant(), false);

        double damping_const = tissue.GetDampingConstant(8);

        TS_ASSERT_DELTA(damping_const, TissueConfig::Instance()->GetDampingConstantNormal(), 1e-6);

        double mutant_damping_const = tissue.GetDampingConstant(9);

        TS_ASSERT_DELTA(mutant_damping_const, TissueConfig::Instance()->GetDampingConstantMutant(), 1e-6);

        tissue.SetAreaBasedDampingConstant(true);

        TS_ASSERT_EQUALS(tissue.UseAreaBasedDampingConstant(), true);

        // Note that this method is usually called by TissueSimulation::Solve()
        tissue.CreateVoronoiTessellation();

        double area_based_damping_const = tissue.GetDampingConstant(8);

        // Since the tissue is in equilibrium, we should get the same damping constant as before
        TS_ASSERT_DELTA(area_based_damping_const, TissueConfig::Instance()->GetDampingConstantNormal(), 1e-6);
    }

    void TestSetNodeAndAddCell()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a tissue, with no ghost nodes at the moment
        MeshBasedTissue<2> tissue(mesh, cells);

        //////////////////
        // Test set node
        //////////////////

        // Move node 0 by a small amount
        AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
        c_vector<double,2> new_location = tissue.GetLocationOfCellCentre(*cell_iter);
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        ChastePoint<2> new_location_point(new_location);
        tissue.SetNode(tissue.GetLocationIndexUsingCell(*cell_iter), new_location_point);

        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[0], new_location[0], 1e-12);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[1], new_location[1], 1e-12);

        //////////////////
        // Test add cell
        //////////////////
        unsigned old_num_nodes = mesh.GetNumNodes();
        unsigned old_num_cells = tissue.rGetCells().size();

        // Create a new cell, DON'T set the node index, set birth time=-1
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        TissueCell cell(STEM, p_state, new FixedDurationGenerationBasedCellCycleModel());
        cell.SetBirthTime(-1);
        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2;
        new_cell_location[1] = 2;

        tissue.AddCell(cell, new_cell_location);

        // Tissue should have updated mesh and cells
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), old_num_nodes+1);

        // Same test via tissue class
        TS_ASSERT_EQUALS(tissue.rGetMesh().GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), old_num_cells+1);

        // Check the location of the new node
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[1], 2.0, 1e-12);

        // Check the index of the new cell
        TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(tissue.rGetCells().back()), old_num_nodes);
    }

    void TestRemoveDeadCellsAndUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        cells[27].StartApoptosis();

        // Create a tissue without ghost nodes
        MeshBasedTissue<2> tissue(mesh,cells);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 81u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 81u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = tissue.RemoveDeadCells();

        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 80u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 80u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 80u);
        TS_ASSERT_DIFFERS(tissue.rGetCells().size(), cells.size()); // Tissue now copies cells

        tissue.Update();

        // For coverage
        NodeMap map(mesh.GetNumAllNodes());
        map.ResetToIdentity();
        tissue.UpdateGhostNodesAfterReMesh(map);

        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 80u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh.GetNumAllNodes());

        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
        {
            TS_ASSERT_EQUALS(tissue.IsGhostNode(i), false);
        }

        // Finally, check the cells node indices have updated

        // We expect the cell node indices to be {0,11,...,79}
        std::set<unsigned> expected_node_indices;
        for (unsigned i=0; i<tissue.GetNumRealCells(); i++)
        {
            expected_node_indices.insert(i);
        }

        // Get actual cell node indices
        std::set<unsigned> node_indices;

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Record node index corresponding to cell
            unsigned node_index = tissue.GetLocationIndexUsingCell(*cell_iter);
            node_indices.insert(node_index);
        }

        TS_ASSERT_EQUALS(node_indices, expected_node_indices);
    }

    void TestUpdateNodeLocations()
    {
        // Test MeshBasedTissue::UpdateNodeLocations()

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a tissue, with no ghost nodes at the moment
        MeshBasedTissue<2> tissue(mesh, cells);

        // Make up some forces
        std::vector<c_vector<double, 2> > old_posns(tissue.GetNumNodes());
        std::vector<c_vector<double, 2> > forces_on_nodes(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            old_posns[i][0] = tissue.GetNode(i)->rGetLocation()[0];
            old_posns[i][1] = tissue.GetNode(i)->rGetLocation()[1];

            forces_on_nodes[i][0] = i*0.01;
            forces_on_nodes[i][1] = 2*i*0.01;
        }

        // Call method
        double time_step = 0.01;
        tissue.UpdateNodeLocations(forces_on_nodes, time_step);

        // Check that node locations were correctly updated
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(tissue.GetNode(i)->rGetLocation()[0], old_posns[i][0] +   i*0.01*0.01, 1e-9);
            TS_ASSERT_DELTA(tissue.GetNode(i)->rGetLocation()[1], old_posns[i][1] + 2*i*0.01*0.01, 1e-9);
        }
    }

    void TestOutputWriters()
    {
        // Resetting the Maximum cell Id to zero (to account for previous tests)
        TissueCell::ResetMaxCellId();

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        cells[0].SetCellProliferativeType(APOPTOTIC); // coverage

        // Cover mutation state reporting
        boost::shared_ptr<AbstractCellMutationState> p_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellMutationState> p_labelled(CellMutationStateRegistry::Instance()->Get<LabelledCellMutationState>());
        boost::shared_ptr<AbstractCellMutationState> p_apc1(CellMutationStateRegistry::Instance()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellMutationState> p_apc2(CellMutationStateRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(CellMutationStateRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>());

        cells[0].SetMutationState(p_state);
        cells[1].SetMutationState(p_apc1);
        cells[2].SetMutationState(p_apc2);
        cells[3].SetMutationState(p_bcat1);
        cells[4].SetMutationState(p_labelled);

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Test set methods

        TissueConfig::Instance()->SetOutputVoronoiData(true);
        TissueConfig::Instance()->SetOutputTissueAreas(true);
        TissueConfig::Instance()->SetOutputCellAreas(true);

        // This method is usually called by Update()
        tissue.CreateVoronoiTessellation();

        std::string output_directory = "TestTissueWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        TissueConfig::Instance()->SetOutputCellMutationStates(true);
        TissueConfig::Instance()->SetOutputCellProliferativeTypes(true);
        TissueConfig::Instance()->SetOutputCellAges(true);

        tissue.CreateOutputFiles(output_directory, false);
        tissue.WriteResultsToFiles();
        tissue.CloseOutputFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizelements   cell_based/test/data/TestTissueWriters/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes      cell_based/test/data/TestTissueWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes  cell_based/test/data/TestTissueWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "tissueareas.dat       cell_based/test/data/TestTissueWriters/tissueareas.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat         cell_based/test/data/TestTissueWriters/cellareas.dat").c_str()), 0);

        // Test the GetCellMutationStateCount function: there should only be healthy cells
        std::vector<unsigned> cell_mutation_states = tissue.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 5u);
        for (unsigned i=0; i<cell_mutation_states.size(); i++)
        {   // There should be one of each kind of mutation in the results files.
            TS_ASSERT_EQUALS(cell_mutation_states[i], 1u);
        }

        // Test the GetCellProliferativeTypeCount function - we should have 4 stem cells and 1 dead cell (for coverage)
        std::vector<unsigned> cell_types = tissue.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 4u);
        TS_ASSERT_EQUALS(cell_types[0], 4u);
        TS_ASSERT_EQUALS(cell_types[1], 0u);
        TS_ASSERT_EQUALS(cell_types[2], 0u);
        TS_ASSERT_EQUALS(cell_types[3], 1u);
    }

    void TestGetLocationOfCellCentreAndGetNodeCorrespondingToCell() throw (Exception)
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

        // Create the tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Loop over nodes
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Record node location
            c_vector<double,2> node_location = tissue.GetLocationOfCellCentre(*cell_iter);

            // Test GetLocationOfCellCentre()
            TS_ASSERT_DELTA(node_location[0], tissue.GetLocationOfCellCentre(*cell_iter)[0], 1e-9);
            TS_ASSERT_DELTA(node_location[1], tissue.GetLocationOfCellCentre(*cell_iter)[1], 1e-9);
        }
    }

    // This test checks that the cells and nodes are correctly archived.
    void TestArchivingMeshBasedTissue() throw (Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mesh_based_tissue.arch";
        ArchiveLocationInfo::SetMeshFilename("mesh_based_tissue_mesh");

        std::vector<c_vector<double,2> > cell_locations;

        // Archive a tissue
        {
            // Need to set up time
            unsigned num_steps=10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // Set up cells, one for each node. Give each a birth time of -node_index,
            // so the age = node_index
            std::vector<TissueCell> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

            // Create the tissue
            MeshBasedTissue<2>* const p_tissue = new MeshBasedTissue<2>(mesh, cells);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            for (AbstractTissue<2>::Iterator cell_iter=p_tissue->Begin();
                 cell_iter!=p_tissue->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
                cell_locations.push_back(p_tissue->GetLocationOfCellCentre(*cell_iter));
            }

            std::set<TissueCell*> cell_pair_0_1 = p_tissue->CreateCellPair(p_tissue->rGetCellUsingLocationIndex(0), p_tissue->rGetCellUsingLocationIndex(1));
            p_tissue->MarkSpring(cell_pair_0_1);

            // Set area-based viscosity
            p_tissue->SetAreaBasedDampingConstant(true);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the tissue to the archive
            (*p_arch) << static_cast<const SimulationTime&> (*p_simulation_time);
            (*p_arch) << p_tissue;
            SimulationTime::Destroy();
            delete p_tissue;
        }

        // Restore tissue
        {
            // Need to set up time
            unsigned num_steps=10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            MeshBasedTissue<2>* p_tissue;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> *p_simulation_time;

            (*p_arch) >> p_tissue;

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0;
            for (AbstractTissue<2>::Iterator cell_iter=p_tissue->Begin();
                 cell_iter!=p_tissue->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(),(double)(counter),1e-7);
                TS_ASSERT_DELTA(p_tissue->GetLocationOfCellCentre(*cell_iter)[0],cell_locations[counter][0],1e-9);
                TS_ASSERT_DELTA(p_tissue->GetLocationOfCellCentre(*cell_iter)[1],cell_locations[counter][1],1e-9);
                counter++;
            }

            // Check the marked spring
            std::set<TissueCell*> cell_pair_0_1 = p_tissue->CreateCellPair(p_tissue->rGetCellUsingLocationIndex(0), p_tissue->rGetCellUsingLocationIndex(1));
            TS_ASSERT_EQUALS(p_tissue->IsMarkedSpring(cell_pair_0_1), true);

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the tissue has been restored
            TS_ASSERT_EQUALS(p_tissue->rGetCells().size(), 5u);

            // Check area-based viscosity is still true
            TS_ASSERT_EQUALS(p_tissue->UseAreaBasedDampingConstant(), true);

            TS_ASSERT_EQUALS(p_tissue->rGetMesh().GetNumNodes(), 5u);

            delete p_tissue;
        }
    }

    void TestSpringMarking()
    {
        // Create a small tissue
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0.5));
        nodes.push_back(new Node<2>(1, false, 1, 0));
        nodes.push_back(new Node<2>(2, false, 1, 1));
        nodes.push_back(new Node<2>(3, false, 2, 0.5));
        nodes.push_back(new Node<2>(4, false, 2, 1.5));

        MutableMesh<2,2> mesh(nodes);

        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedTissue<2> tissue(mesh, cells);

        std::set<TissueCell*> cell_pair_1_2 = tissue.CreateCellPair(tissue.rGetCellUsingLocationIndex(1), tissue.rGetCellUsingLocationIndex(2));
        std::set<TissueCell*> cell_pair_3_4 = tissue.CreateCellPair(tissue.rGetCellUsingLocationIndex(3), tissue.rGetCellUsingLocationIndex(4));
        std::set<TissueCell*> cell_pair_1_4 = tissue.CreateCellPair(tissue.rGetCellUsingLocationIndex(1), tissue.rGetCellUsingLocationIndex(4));
        std::set<TissueCell*> cell_pair_0_2 = tissue.CreateCellPair(tissue.rGetCellUsingLocationIndex(0), tissue.rGetCellUsingLocationIndex(2));

        // Mark some springs
        tissue.MarkSpring(cell_pair_1_2);

        // Unmark and re-mark spring (for coverage)
        tissue.UnmarkSpring(cell_pair_1_2);
        tissue.MarkSpring(cell_pair_1_2);

        tissue.MarkSpring(cell_pair_3_4);

        // Check if springs are marked
        TS_ASSERT_EQUALS(tissue.IsMarkedSpring(cell_pair_1_2), true);
        TS_ASSERT_EQUALS(tissue.IsMarkedSpring(cell_pair_3_4), true);

        TS_ASSERT_EQUALS(tissue.IsMarkedSpring(cell_pair_1_4), false);
        TS_ASSERT_EQUALS(tissue.IsMarkedSpring(cell_pair_0_2), false);

        // Delete cell 4
        tissue.rGetCellUsingLocationIndex(4).Kill();
        tissue.RemoveDeadCells();

        // Check springs with non-deleted cells are still marked
        TS_ASSERT_EQUALS(tissue.IsMarkedSpring(cell_pair_1_2), true);
        tissue.CheckTissueCellPointers();

        // Move cell 2
        AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
        ++cell_iter;
        ++cell_iter;
        TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(*cell_iter), 2u);
        ChastePoint<2> new_location(1, 10);
        tissue.SetNode(tissue.GetLocationIndexUsingCell(*cell_iter), new_location);

        // Update tissue
        tissue.Update();

        tissue.CheckTissueCellPointers();

        // Check there is no marked spring between nodes 1 & 2
        TS_ASSERT_EQUALS(tissue.IsMarkedSpring(cell_pair_1_2), false);
    }

    void TestSettingCellAncestors() throw (Exception)
    {
        // Create a small tissue
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0.5));
        nodes.push_back(new Node<2>(1, false, 1, 0));
        nodes.push_back(new Node<2>(2, false, 1, 1));
        nodes.push_back(new Node<2>(3, false, 2, 0.5));
        nodes.push_back(new Node<2>(4, false, 2, 1.5));

        MutableMesh<2,2> mesh(nodes);

        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedTissue<2> tissue(mesh, cells);

        // Test that the tissue makes all cells fix the node index as ancestor
        tissue.SetCellAncestorsToLocationIndices();

        unsigned counter = 0;
        for (AbstractTissue<2>::Iterator cell_iter=tissue.Begin();
             cell_iter!=tissue.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetAncestor(), tissue.GetLocationIndexUsingCell(*cell_iter));
            counter ++;
        }
        TS_ASSERT_EQUALS(counter, 5u);

        // Test that we can recover remaining number of ancestors
        std::set<unsigned> remaining_ancestors = tissue.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 5u);

        // Test that the set correctly represents a monoclonal population
        for (AbstractTissue<2>::Iterator cell_iter=tissue.Begin();
             cell_iter!=tissue.End();
             ++cell_iter)
        {
            // Set all cells to have the same ancestor
            cell_iter->SetAncestor(1u);
        }
        remaining_ancestors = tissue.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 1u);
    }

    void TestIsCellAssociatedWithADeletedLocation() throw (Exception)
    {
        // Create a simple mesh
        HoneycombMeshGenerator generator(4, 4, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->GetNode(0)->MarkAsDeleted();

        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Create tissue but do not try to validate
        MeshBasedTissue<2> tissue(*p_mesh, cells, std::vector<unsigned>(), false, false);

        // Test IsCellAssociatedWithADeletedLocation() method
        for (MeshBasedTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            bool is_deleted = tissue.IsCellAssociatedWithADeletedLocation(*cell_iter);

            if (tissue.GetLocationIndexUsingCell(*cell_iter) == 0)
            {
                TS_ASSERT_EQUALS(is_deleted, true);
            }
            else
            {
                TS_ASSERT_EQUALS(is_deleted, false);
            }
        }
    }
};


#endif /*TESTMESHBASEDTISSUE_HPP_*/
