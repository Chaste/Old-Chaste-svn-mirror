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

#include "CellwiseData.hpp"
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
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"

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
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, DIM> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create the tissue
        unsigned num_cells = cells.size();
        MeshBasedTissue<DIM> tissue(mesh, cells);

        TS_ASSERT_EQUALS(tissue.rGetMesh().GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), num_cells);

        // Test set/get method of member variable
        TS_ASSERT_DELTA(tissue.GetMeinekeDivisionSeparation(), 0.3, 1e-6);
        tissue.SetMeinekeDivisionSeparation(0.5);
        TS_ASSERT_DELTA(tissue.GetMeinekeDivisionSeparation(), 0.5, 1e-6);
        tissue.SetMeinekeDivisionSeparation(0.3);

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
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<mesh.GetNumNodes()-1; i++)
        {
            AbstractCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
            p_cell_cycle_model->SetCellProliferativeType(STEM);
            TissueCellPtr p_cell(new TissueCell(p_state, p_cell_cycle_model));
            double birth_time = 0.0 - i;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Fails as no cell corresponding to node 4
        std::vector<TissueCellPtr> cells_copy(cells);
        TS_ASSERT_THROWS_THIS(MeshBasedTissue<2> tissue2(mesh, cells_copy),
                              "Node 4 does not appear to have a cell associated with it");

        // Add another cell
        AbstractCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_cycle_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_cell(new TissueCell(p_state, p_cell_cycle_model));
        double birth_time = -4.0;
        p_cell->SetBirthTime(birth_time);
        cells.push_back(p_cell);

        std::vector<TissueCellPtr> cells_copy2(cells);
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
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Give cells 0 and 1 specific mutations to enable later testing
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        boost::shared_ptr<AbstractCellProperty> p_apc1(CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>());
        cells[0]->AddCellProperty(p_label);
        cells[1]->SetMutationState(p_apc1);

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Create cell pair
        std::set<TissueCellPtr> cell_pair = tissue.CreateCellPair(cells[0], cells[1]);

        // Check the cell pair was created correctly
        std::set<TissueCellPtr>::iterator cell_pair_iter = cell_pair.begin();

        TissueCellPtr p_cell0 = *cell_pair_iter;
        TS_ASSERT_EQUALS(p_cell0->HasCellProperty<CellLabel>(), true);

        ++cell_pair_iter;
        TissueCellPtr p_cell1 = *cell_pair_iter;
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
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Bestow mutations on some cells
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_apc1(new ApcOneHitCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_apc2(new ApcTwoHitCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_bcat1(new BetaCateninOneHitCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);

        cells[0]->SetMutationState(p_state);
        cells[1]->SetMutationState(p_apc1);
        cells[2]->SetMutationState(p_apc2);
        cells[3]->SetMutationState(p_bcat1);
        cells[4]->AddCellProperty(p_label);

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

        // Coverage
        TS_ASSERT_DELTA(tissue.GetAreaBasedDampingConstantParameter(), 0.1, 1e-6);
        tissue.SetAreaBasedDampingConstantParameter(0.5);
        TS_ASSERT_DELTA(tissue.GetAreaBasedDampingConstantParameter(), 0.5, 1e-6);

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
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        boost::shared_ptr<AbstractCellProperty> p_apc2(new ApcTwoHitCellMutationState);
        cells[9]->SetMutationState(p_apc2);

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

        std::vector<TissueCellPtr> cells;
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
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);

        FixedDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_cycle_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_cell(new TissueCell(p_state, p_cell_cycle_model));
        p_cell->SetBirthTime(-1);
        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2;
        new_cell_location[1] = 2;

        tissue.AddCell(p_cell, new_cell_location);

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
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        cells[27]->StartApoptosis();

        // Create a tissue without ghost nodes
        MeshBasedTissue<2> tissue(mesh, cells);

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

        std::vector<TissueCellPtr> cells;
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

    void TestVoronoiMethods()
    {
        // First test the 2D case...

        // Create 2D mesh
        std::vector<Node<2> *> nodes2d;
        nodes2d.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes2d.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes2d.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes2d.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes2d.push_back(new Node<2>(4, false, 0.5, 0.5));
        MutableMesh<2,2> mesh2d(nodes2d);

        // Create cells
        std::vector<TissueCellPtr> cells2d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator2d;
        cells_generator2d.GenerateBasic(cells2d, mesh2d.GetNumNodes());

        // Create tissue
        MeshBasedTissue<2> tissue2d(mesh2d, cells2d);

        // Create Voronoi tessellation
        tissue2d.CreateVoronoiTessellation();

        // Test element areas
        TS_ASSERT_DELTA(tissue2d.GetVolumeOfVoronoiElement(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(tissue2d.GetVolumeOfVoronoiElement(1), 0.0, 1e-6);
        TS_ASSERT_DELTA(tissue2d.GetVolumeOfVoronoiElement(2), 0.0, 1e-6);
        TS_ASSERT_DELTA(tissue2d.GetVolumeOfVoronoiElement(3), 0.0, 1e-6);
        TS_ASSERT_DELTA(tissue2d.GetVolumeOfVoronoiElement(4), 0.5, 1e-6);

        // Test element perimeters
        TS_ASSERT_DELTA(tissue2d.GetSurfaceAreaOfVoronoiElement(0), sqrt(2), 1e-6);
        TS_ASSERT_DELTA(tissue2d.GetSurfaceAreaOfVoronoiElement(1), sqrt(2), 1e-6);
        TS_ASSERT_DELTA(tissue2d.GetSurfaceAreaOfVoronoiElement(2), sqrt(2), 1e-6);
        TS_ASSERT_DELTA(tissue2d.GetSurfaceAreaOfVoronoiElement(3), sqrt(2), 1e-6);
        TS_ASSERT_DELTA(tissue2d.GetSurfaceAreaOfVoronoiElement(4), 2*sqrt(2), 1e-6);

        // ...now test the 3D case

        // Create 3D mesh
        std::vector<Node<3>*> nodes3d;
        nodes3d.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes3d.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes3d.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes3d.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh3d(nodes3d);

        // Create cells
        std::vector<TissueCellPtr> cells3d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator3d;
        cells_generator3d.GenerateBasic(cells3d, mesh3d.GetNumNodes());

        // Create tissue
        MeshBasedTissue<3> tissue3d(mesh3d, cells3d);

        // Create Voronoi tessellation
        tissue3d.CreateVoronoiTessellation();

        // Test element volumes and surface areas
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_THROWS_THIS(tissue3d.GetVolumeOfVoronoiElement(i),
                                  "This index does not correspond to a VertexElement");

            TS_ASSERT_THROWS_THIS(tissue3d.GetSurfaceAreaOfVoronoiElement(i),
                                  "This index does not correspond to a VertexElement");
        }

        // The Voronoi tessellation should comprise a single tetrahedral VertexElement
        TS_ASSERT_EQUALS(tissue3d.GetVoronoiTessellation()->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(tissue3d.GetVoronoiTessellation()->GetNumFaces(), 4u);
        TS_ASSERT_EQUALS(tissue3d.GetVoronoiTessellation()->GetNumElements(), 1u);

        // The faces are not all equal
        for (unsigned face_index=0; face_index<4; face_index++)
        {
            VertexElement<2,3>* p_face = tissue3d.GetVoronoiTessellation()->GetFace(face_index);

            if (face_index == 1)
            {
                TS_ASSERT_DELTA(tissue3d.GetVoronoiTessellation()->GetAreaOfFace(p_face), 1.9485, 1e-4);
            }
            else
            {
                TS_ASSERT_DELTA(tissue3d.GetVoronoiTessellation()->GetAreaOfFace(p_face), 1.125, 1e-4);
            }
        }

        TS_ASSERT_DELTA(tissue3d.GetVolumeOfVoronoiElement(4), 0.6495, 1e-4);
        TS_ASSERT_DELTA(tissue3d.GetSurfaceAreaOfVoronoiElement(4), 3*1.125 + 1.9485, 1e-4);

        // Check that the Voronoi tessellation can be returned successfully as a reference
        VertexMesh<3,3>* p_tessellation1 = tissue3d.GetVoronoiTessellation();
        TS_ASSERT_EQUALS(p_tessellation1->GetNumNodes(), 4u);

        // Move node 0 by a small amount
        AbstractTissue<3>::Iterator cell_iter = tissue3d.Begin();
        c_vector<double,3> new_location = tissue3d.GetLocationOfCellCentre(*cell_iter);
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        new_location[2] += 1e-2;
        ChastePoint<3> new_location_point(new_location);
        tissue3d.SetNode(tissue3d.GetLocationIndexUsingCell(*cell_iter), new_location_point);

        // Re-create Voronoi tessellation
        tissue3d.CreateVoronoiTessellation();

        VertexMesh<3,3>* p_tessellation2 = tissue3d.GetVoronoiTessellation();
        TS_ASSERT_EQUALS(p_tessellation2->GetNumNodes(), 4u);

//        TS_ASSERT_EQUALS(p_tessellation1->GetNumNodes(), 4u);
    }

    void TestTissueWritersIn2d()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Resetting the Maximum cell Id to zero (to account for previous tests)
        TissueCell::ResetMaxCellId();

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Cover mutation state reporting
        boost::shared_ptr<AbstractCellProperty> p_apc1(CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc2(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(CellPropertyRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        cells[0]->AddCellProperty(p_apoptotic_state);
        cells[1]->SetMutationState(p_apc1);
        cells[2]->SetMutationState(p_apc2);
        cells[3]->SetMutationState(p_bcat1);
        cells[4]->AddCellProperty(p_label);

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        TS_ASSERT_EQUALS(tissue.GetIdentifier(), "MeshBasedTissue<2>");

        // Test set/get methods
        TS_ASSERT_EQUALS(tissue.GetOutputVoronoiData(), false);
        TS_ASSERT_EQUALS(tissue.GetOutputCellIdData(), false);

        tissue.SetOutputVoronoiData(true);
        tissue.SetOutputCellIdData(true);

        TS_ASSERT_EQUALS(tissue.GetOutputVoronoiData(), true);
        TS_ASSERT_EQUALS(tissue.GetOutputCellIdData(), true);

        // Coverage of writing CellwiseData to VTK
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(tissue.GetNumRealCells(), 2);
        p_data->SetTissue(&tissue);
        for (unsigned var=0; var<2; var++)
        {
            for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
                 cell_iter != tissue.End();
                 ++cell_iter)
            {
                p_data->SetValue((double) 3.0*var, tissue.GetLocationIndexUsingCell(*cell_iter), var);
            }
        }

        // Test set methods
        tissue.SetOutputTissueVolumes(true);
        tissue.SetOutputCellVolumes(true);
        tissue.SetOutputCellAncestors(true);
        tissue.SetOutputCellMutationStates(true);
        tissue.SetOutputCellProliferativeTypes(true);
        tissue.SetOutputCellAges(true);
        tissue.SetOutputCellCyclePhases(true);

        // This method is usually called by Update()
        tissue.CreateVoronoiTessellation();

        std::string output_directory = "TestTissueWritersIn2d";
        OutputFileHandler output_file_handler(output_directory, false);

        tissue.CreateOutputFiles(output_directory, false);
        tissue.WriteResultsToFiles();
        tissue.CloseOutputFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizelements   cell_based/test/data/TestTissueWritersIn2d/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes      cell_based/test/data/TestTissueWritersIn2d/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes  cell_based/test/data/TestTissueWritersIn2d/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "tissueareas.dat       cell_based/test/data/TestTissueWritersIn2d/tissueareas.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat         cell_based/test/data/TestTissueWritersIn2d/cellareas.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat         cell_based/test/data/TestTissueWritersIn2d/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "voronoi.dat           cell_based/test/data/TestTissueWritersIn2d/voronoi.dat").c_str()), 0);

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = tissue.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 2u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 1u);

        // Test the GetCellProliferativeTypeCount function - we should have 4 stem cells and 1 dead cell (for coverage)
        std::vector<unsigned> cell_types = tissue.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 3u);
        TS_ASSERT_EQUALS(cell_types[0], 5u);
        TS_ASSERT_EQUALS(cell_types[1], 0u);
        TS_ASSERT_EQUALS(cell_types[2], 0u);

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestTissueWritersIn3d()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Resetting the Maximum cell Id to zero (to account for previous tests)
        TissueCell::ResetMaxCellId();

        // Create a simple 3D mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        cells[4]->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>()); // coverage

        // Create tissue
        MeshBasedTissue<3> tissue(mesh, cells);

        TS_ASSERT_EQUALS(tissue.GetIdentifier(), "MeshBasedTissue<3>");

        // Test set methods
        tissue.SetOutputVoronoiData(true);
        tissue.SetOutputTissueVolumes(true);
        tissue.SetOutputCellVolumes(true);
        tissue.SetOutputCellAncestors(true);
        tissue.SetOutputCellMutationStates(true);
        tissue.SetOutputCellProliferativeTypes(true);
        tissue.SetOutputCellAges(true);
        tissue.SetOutputCellCyclePhases(true);

        // This method is usually called by Update()
        tissue.CreateVoronoiTessellation();

        std::string output_directory = "TestTissueWritersIn3d";
        OutputFileHandler output_file_handler(output_directory, false);

        tissue.CreateOutputFiles(output_directory, false);
        tissue.WriteResultsToFiles();
        tissue.CloseOutputFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizelements   cell_based/test/data/TestTissueWritersIn3d/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes      cell_based/test/data/TestTissueWritersIn3d/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes  cell_based/test/data/TestTissueWritersIn3d/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "tissueareas.dat       cell_based/test/data/TestTissueWritersIn3d/tissueareas.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat         cell_based/test/data/TestTissueWritersIn3d/cellareas.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat         cell_based/test/data/TestTissueWritersIn3d/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "voronoi.dat           cell_based/test/data/TestTissueWritersIn3d/voronoi.dat").c_str()), 0);

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = tissue.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 5u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 0u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 0u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 0u);

        // Test the GetCellProliferativeTypeCount function
        std::vector<unsigned> cell_types = tissue.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 3u);
        TS_ASSERT_EQUALS(cell_types[0], 5u);
        TS_ASSERT_EQUALS(cell_types[1], 0u);
        TS_ASSERT_EQUALS(cell_types[2], 0u);
    }

    void TestGetLocationOfCellCentreAndGetNodeCorrespondingToCell() throw (Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCellPtr> cells;
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
            std::vector<TissueCellPtr> cells;
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

            std::set<TissueCellPtr> cell_pair_0_1 = p_tissue->CreateCellPair(p_tissue->GetCellUsingLocationIndex(0), p_tissue->GetCellUsingLocationIndex(1));
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
                TS_ASSERT_DELTA(p_tissue->GetLocationOfCellCentre(*cell_iter)[0], cell_locations[counter][0], 1e-9);
                TS_ASSERT_DELTA(p_tissue->GetLocationOfCellCentre(*cell_iter)[1], cell_locations[counter][1], 1e-9);
                counter++;
            }

            // Check the marked spring
            std::set<TissueCellPtr> cell_pair_0_1 = p_tissue->CreateCellPair(p_tissue->GetCellUsingLocationIndex(0), p_tissue->GetCellUsingLocationIndex(1));
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

        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedTissue<2> tissue(mesh, cells);

        std::set<TissueCellPtr> cell_pair_1_2 = tissue.CreateCellPair(tissue.GetCellUsingLocationIndex(1), tissue.GetCellUsingLocationIndex(2));
        std::set<TissueCellPtr> cell_pair_3_4 = tissue.CreateCellPair(tissue.GetCellUsingLocationIndex(3), tissue.GetCellUsingLocationIndex(4));
        std::set<TissueCellPtr> cell_pair_1_4 = tissue.CreateCellPair(tissue.GetCellUsingLocationIndex(1), tissue.GetCellUsingLocationIndex(4));
        std::set<TissueCellPtr> cell_pair_0_2 = tissue.CreateCellPair(tissue.GetCellUsingLocationIndex(0), tissue.GetCellUsingLocationIndex(2));

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
        tissue.GetCellUsingLocationIndex(4)->Kill();
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

        std::vector<TissueCellPtr> cells;
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
        std::vector<TissueCellPtr> cells;
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
