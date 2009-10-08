/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTFORCES_HPP_
#define TESTFORCES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "FixedDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "IngeWntSwatCellCycleModelCellsGenerator.hpp"
#include "MeshBasedTissueWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "WntConcentration.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestForces : public AbstractCancerTestSuite
{
public:

    void TestGeneralisedLinearSpringForceMethods() throw (Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<TissueCell> cells;
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            CryptCellMutationState mutation_state = HEALTHY;
            if (i==60)
            {
                mutation_state = APC_TWO_HIT;
            }

            TissueCell cell(STEM, mutation_state, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(-10);
            cells.push_back(cell);
        }

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);
        GeneralisedLinearSpringForce<2> linear_force;

        /*
         ************************************************************************
         ************************************************************************
         *  Test node force calculation
         ************************************************************************
         ************************************************************************
         */

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }

        linear_force.AddForceContribution(node_forces, tissue);

        // Test forces on non-ghost nodes
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            unsigned node_index = tissue.GetLocationIndexUsingCell(*cell_iter);

            TS_ASSERT_DELTA(node_forces[node_index][0], 0.0, 1e-4);
            TS_ASSERT_DELTA(node_forces[node_index][1], 0.0, 1e-4);
        }

        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point = p_mesh->GetNode(59)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];

        p_mesh->SetNode(59, new_point, false);

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > new_node_forces;
        new_node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             new_node_forces.push_back(zero_vector<double>(2));
        }
        linear_force.AddForceContribution(new_node_forces, tissue);

        TS_ASSERT_DELTA(new_node_forces[60][0], 0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[60][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[59][0], (-3+4.0/sqrt(7))*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[59][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[58][0], 0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[58][1], 0.0, 1e-4);

        /*
         ************************************************************************
         ************************************************************************
         *  Test spring force calculation
         ************************************************************************
         ************************************************************************
         */

        c_vector<double,2> force_on_spring; // between nodes 59 and 60

        // Find one of the elements that nodes 59 and 60 live on
        ChastePoint<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0] + 0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01;

        unsigned elem_index = p_mesh->GetContainingElementIndex(new_point2, false);
        Element<2,2>* p_element = p_mesh->GetElement(elem_index);

        force_on_spring = linear_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                   p_element->GetNodeGlobalIndex(0),
                                                                   tissue);

        TS_ASSERT_DELTA(force_on_spring[0], 0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);

        /*
         ******************************************
         ******************************************
         *  Test force with cutoff point
         ******************************************
         ******************************************
         */
        double dist = norm_2( p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(),
                              p_element->GetNode(1)->rGetLocation()) );

        linear_force.UseCutoffPoint(dist-0.1);

        // Coverage
        TS_ASSERT_DELTA(linear_force.GetCutoffPoint(), dist-0.1, 1e-4);

        force_on_spring = linear_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                   p_element->GetNodeGlobalIndex(0),
                                                                   tissue);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);
    }


    void TestGeneralisedLinearSpringForceWithEdgeLengthBasedSpring() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false, crypt_width/cells_across);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true); // true = mature cells

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);
        LinearSpringWithVariableSpringConstantsForce<2> linear_force;

        // Check that the force between nodes is correctly calculated when the 'spring constant' is constant
        linear_force.SetEdgeBasedSpringConstant(false);

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
            spring_iterator != tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                 nodeB_global_index,
                                                                                 tissue);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 6.25, 1e-3);
        }

        // Check that the force between nodes is correctly calculated when the 'spring constant'
        // is proportional to the length of the edge between adjacent cells
        linear_force.SetEdgeBasedSpringConstant(true);
        tissue.CreateVoronoiTessellation();  // normally done in a simulation loop

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
             spring_iterator != tissue.SpringsEnd();
             ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                 nodeB_global_index,
                                                                                 tissue);

            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 4.34027778, 1e-3);
        }

        // Choose two interior neighbour nodes
        c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(41, 42, tissue);
        TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 4.34027778, 1e-3);

        // Now move node 42 a bit and check that the force calculation changes correctly
        c_vector<double,2> shift;
        shift[0] = 0.1;
        shift[1] = 0.0;
        ChastePoint<2> new_point(p_mesh->GetNode(42u)->rGetLocation() + shift);
        p_mesh->SetNode(21, new_point, false);

        // Check that the new force between nodes is correctly calculated
        tissue.CreateVoronoiTessellation();
        c_vector<double, 2> new_force = linear_force.CalculateForceBetweenNodes(41, 42, tissue);

        // Force calculation: shift is along x-axis so we should have
        // new_edge_length = (5/6 + shift[0])*tan(0.5*arctan(5*sqrt(3)/(5 + 12*shift[0]))),
        // force^2 = mu^2 * (new_edge_length*sqrt(3))^2 * (1 - 5/6 - shift[0])^2
        TS_ASSERT_DELTA(new_force[0]*new_force[0] + new_force[1]*new_force[1], 4.34024, 1e-3);
    }

    // Test on a periodic mesh
    void TestGeneralisedLinearSpringForceWithEdgeBasedSpringsOnPeriodicMesh() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);// true = mature cells

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);
        LinearSpringWithVariableSpringConstantsForce<2> linear_force;

        // Check that the force between nodes is correctly calculated when the spring constant is constant (!)
        linear_force.SetEdgeBasedSpringConstant(false);

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
             spring_iterator != tissue.SpringsEnd();
             ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                 nodeB_global_index,
                                                                                 tissue);

            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 6.25, 1e-3);
        }

        // Check that the force between nodes is correctly calculated when the spring constant
        // is proportional to the length of the edge between adjacenet cells
        linear_force.SetEdgeBasedSpringConstant(true);
        tissue.CreateVoronoiTessellation();

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
             spring_iterator != tissue.SpringsEnd();
             ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                 nodeB_global_index,
                                                                                 tissue);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 4.34027778, 1e-3);
        }
    }

    void TestGeneralisedLinearSpringForceWithSpringConstantsForMutantCells()
    {
        // Create a small tissue
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0));
        nodes.push_back(new Node<2>(1, false, 0, 2));
        nodes.push_back(new Node<2>(2, false, 2, 2));
        nodes.push_back(new Node<2>(3, false, 2, 0));

        MutableMesh<2,2> mesh(nodes);

        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedTissue<2> tissue(mesh, cells);

        LinearSpringWithVariableSpringConstantsForce<2> linear_force;

        // Set cells mutation states
        tissue.rGetCellUsingLocationIndex(0).SetMutationState(HEALTHY);
        tissue.rGetCellUsingLocationIndex(1).SetMutationState(LABELLED);
        tissue.rGetCellUsingLocationIndex(2).SetMutationState(APC_TWO_HIT);
        tissue.rGetCellUsingLocationIndex(3).SetMutationState(BETA_CATENIN_ONE_HIT);

        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(0, 1, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(1, 2, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(2, 3, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(3, 0, tissue)), 15.0, 1e-10);

        linear_force.SetMutantSprings(true);

        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(0, 1, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(1, 2, tissue)), 22.5, 1e-10);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(2, 3, tissue)), 30.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(3, 0, tissue)), 22.5, 1e-10);

        linear_force.SetMutantSprings(true, 4.0, 3.0);

        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(0, 1, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(1, 2, tissue)), 45.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(2, 3, tissue)), 60.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(3, 0, tissue)), 45.0, 1e-10);
    }


    void TestGeneralisedLinearSpringForceWithSpringConstantsForIngeBCatCells()
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(6, 12, 0, true, 1.1);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        IngeWntSwatCellCycleModelCellsGenerator<2> cells_generator(2u);
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, false);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // As there is no tissue simulation, we must explicitly initialise the cells
        crypt.InitialiseCells();

        LinearSpringWithVariableSpringConstantsForce<2> linear_force;

        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.50, 1e-10);

        linear_force.SetBetaCateninSprings(true);
        crypt.CreateVoronoiTessellation();  // this method is normally called in a simulation loop

        /// \todo this is currently a rather poor test - it just checks that
        /// there is SOME dependency of the spring constant on the beta catenin level
        /// experienced by both cells (see #627) - we need to get rid of the magic numbers!

        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.5*8.59312/18.14, 1e-5);

        TissueConfig::Instance()->SetBetaCatSpringScaler(20/6.0);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.5*8.59312/20.0, 1e-5);

        WntConcentration<2>::Destroy();
    }

    void TestGeneralisedLinearSpringForceWithSpringConstantsForApoptoticCells()
    {
        // Set up stretched tissue
        HoneycombMeshGenerator generator(4, 4, 0, false, 2.0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        MeshBasedTissueWithGhostNodes<2> stretched_tissue(*p_mesh, cells, location_indices);

        // As there is no tissue simulation we must explicitly initialise the cells
        stretched_tissue.InitialiseCells();

        // Set one of the non-boundary cells to be necrotic
        stretched_tissue.rGetCellUsingLocationIndex(6).SetCellProliferativeType(APOPTOTIC);

        LinearSpringWithVariableSpringConstantsForce<2> linear_force;
        linear_force.SetApoptoticSprings(true);

        TS_ASSERT_EQUALS( stretched_tissue.rGetCellUsingLocationIndex(6).GetCellProliferativeType(), APOPTOTIC);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(6, 10, stretched_tissue)), 3.3333, 1e-4);

        // Set a neighbouring cell to be necrotic
        stretched_tissue.rGetCellUsingLocationIndex(10).SetCellProliferativeType(APOPTOTIC);

        TS_ASSERT_EQUALS( stretched_tissue.rGetCellUsingLocationIndex(10).GetCellProliferativeType(), APOPTOTIC);
        TS_ASSERT_DELTA( norm_2(linear_force.CalculateForceBetweenNodes(6, 10, stretched_tissue)), 1.8750, 1e-4);

        // Now do similar tests for a squashed tissue
        HoneycombMeshGenerator generator2(4, 4, 0, false, 0.5);
        MutableMesh<2,2>* p_mesh2 = generator2.GetMesh();
        std::vector<unsigned> location_indices2 = generator2.GetCellLocationIndices();

        std::vector<TissueCell> cells2;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator2;
        cells_generator2.GenerateGivenLocationIndices(cells2, location_indices2);

        MeshBasedTissueWithGhostNodes<2> squashed_tissue(*p_mesh2, cells2, location_indices2);
        squashed_tissue.InitialiseCells();

        squashed_tissue.rGetCellUsingLocationIndex(6).SetCellProliferativeType(APOPTOTIC);

        LinearSpringWithVariableSpringConstantsForce<2> linear_force2;
        linear_force2.SetApoptoticSprings(true);

        TS_ASSERT_DELTA( norm_2(linear_force2.CalculateForceBetweenNodes(6, 10, squashed_tissue)), 4.0909, 1e-4);

        squashed_tissue.rGetCellUsingLocationIndex(10).SetCellProliferativeType(APOPTOTIC);

        TS_ASSERT_DELTA( norm_2(linear_force2.CalculateForceBetweenNodes(6, 10, squashed_tissue)), 2.8125, 1e-4);
    }


    void TestGeneralisedLinearSpringForceCalculationIn1d() throw (Exception)
    {
        // Create a 1D mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(5);

        // Set up cells
        std::vector<TissueCell> cells;
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - node_index;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        MeshBasedTissue<1> tissue(mesh, cells);

        // Create force law object
        GeneralisedLinearSpringForce<1> linear_force;

        // Initialise a vector of node forces
        std::vector<c_vector<double, 1> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(1));
        }

        // Compute forces on nodes
        linear_force.AddForceContribution(node_forces, tissue);

        // Test that all springs are in equilibrium
        for (unsigned node_index=0; node_index<tissue.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(node_forces[node_index](0), 0.0, 1e-6);
        }

        // Scale entire mesh and check that forces are correctly calculated
        double scale_factor = 1.5;
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            c_vector<double,1> old_point = mesh.GetNode(node_index)->rGetLocation();
            ChastePoint<1> new_point;
            new_point.rGetLocation()[0] = scale_factor*old_point[0];
            mesh.SetNode(node_index, new_point, false);
        }

        // Recalculate node forces (we can re-use node_forces
        // as previously each node had zero net force on it)
        linear_force.AddForceContribution(node_forces, tissue);

        TissueConfig* p_params = TissueConfig::Instance();
        for (unsigned node_index=0; node_index<tissue.GetNumNodes(); node_index++)
        {
            if (node_index == 0)
            {
                // The first node only experiences a force from its neighbour to the right
                TS_ASSERT_DELTA(node_forces[node_index](0), p_params->GetSpringStiffness()*(scale_factor-1), 1e-6);
            }
            else if (node_index == tissue.GetNumNodes()-1)
            {
                // The last node only experiences a force from its neighbour to the left
                TS_ASSERT_DELTA(node_forces[node_index](0), -p_params->GetSpringStiffness()*(scale_factor-1), 1e-6);
            }
            else
            {
                // The net force on each interior node should still be zero
                TS_ASSERT_DELTA(node_forces[node_index](0), 0.0, 1e-6);
            }
        }

        // Create another tissue and force law
        MutableMesh<1,1> mesh2;
        mesh2.ConstructLinearMesh(5);

        MeshBasedTissue<1> tissue2(mesh2, cells);
        GeneralisedLinearSpringForce<1> linear_force2;

        // Move one node and check that forces are correctly calculated
        ChastePoint<1> shifted_point;
        shifted_point.rGetLocation()[0] = 2.5;
        mesh2.SetNode(2, shifted_point);

        c_vector<double,1> force_between_1_and_2 = linear_force2.CalculateForceBetweenNodes(1, 2, tissue2);
        TS_ASSERT_DELTA(force_between_1_and_2[0], p_params->GetSpringStiffness()*0.5, 1e-6);

        c_vector<double,1> force_between_2_and_3 = linear_force2.CalculateForceBetweenNodes(2, 3, tissue2);
        TS_ASSERT_DELTA(force_between_2_and_3[0], -p_params->GetSpringStiffness()*0.5, 1e-6);

        // Initialise a vector of node forces
        std::vector<c_vector<double,1> > node_forces2;
        node_forces2.reserve(tissue2.GetNumNodes());

        for (unsigned i=0; i<tissue2.GetNumNodes(); i++)
        {
             node_forces2.push_back(zero_vector<double>(1));
        }

        linear_force2.AddForceContribution(node_forces2, tissue2);

        TS_ASSERT_DELTA(node_forces2[2](0), -p_params->GetSpringStiffness(), 1e-6);
    }


    void TestGeneralisedLinearSpringForceCalculationIn3d() throw (Exception)
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;
        TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            cell.SetBirthTime(-50.0);
            cells.push_back(cell);
        }

        MeshBasedTissue<3> tissue(mesh,cells);
        GeneralisedLinearSpringForce<3> linear_force;

        // Test forces on springs
        unsigned nodeA = 0, nodeB = 1;
        Element<3,3>* p_element = mesh.GetElement(0);
        c_vector<double, 3> force = linear_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(nodeA),
                                                                            p_element->GetNodeGlobalIndex(nodeB),
                                                                            tissue);
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(force[i], 0.0, 1e-6);
        }

        // Initialise a vector of node forces
        std::vector<c_vector<double, 3> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(3));
        }

        linear_force.AddForceContribution(node_forces, tissue);

        for (unsigned j=0; j<4; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(node_forces[j](k), 0.0, 1e-6);
            }
        }

        // Scale entire mesh and check that forces are correctly calculated
        TissueConfig* p_params = TissueConfig::Instance();
        double scale_factor = 1.5;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> old_point = mesh.GetNode(i)->rGetLocation();
            ChastePoint<3> new_point;
            new_point.rGetLocation()[0] = scale_factor*old_point[0];
            new_point.rGetLocation()[1] = scale_factor*old_point[1];
            new_point.rGetLocation()[2] = scale_factor*old_point[2];
            mesh.SetNode(i, new_point, false);
        }

        // Recalculate node forces (we can just re-use node_forces,
        // as previously each node had zero net force on it)
        linear_force.AddForceContribution(node_forces, tissue);

        for (unsigned j=0; j<4; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(fabs(node_forces[j](k)), p_params->GetSpringStiffness()*(scale_factor-1)*sqrt(2),1e-6);
            }
        }

        // Move one node and check that forces are correctly calculated
        MutableMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);

        MeshBasedTissue<3> tissue2(mesh2,cells);
        GeneralisedLinearSpringForce<3> linear_force2;

        c_vector<double,3> old_point = mesh2.GetNode(0)->rGetLocation();
        ChastePoint<3> new_point;
        new_point.rGetLocation()[0] = 0.0;
        new_point.rGetLocation()[1] = 0.0;
        new_point.rGetLocation()[2] = 0.0;
        mesh2.SetNode(0, new_point, false);

        unsigned nodeA2 = 0, nodeB2 = 1;
        Element<3,3>* p_element2 = mesh2.GetElement(0);
        c_vector<double,3> force2 = linear_force2.CalculateForceBetweenNodes(p_element2->GetNodeGlobalIndex(nodeA2),
                                                                             p_element2->GetNodeGlobalIndex(nodeB2),
                                                                             tissue2);

        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(fabs(force2[i]),p_params->GetSpringStiffness()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }

        // Initialise a vector of node forces
        std::vector<c_vector<double,3> > node_forces2;
        node_forces2.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces2.push_back(zero_vector<double>(3));
        }

        linear_force2.AddForceContribution(node_forces2, tissue2);

        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(node_forces2[0](i),p_params->GetSpringStiffness()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
    }


    void TestGeneralisedLinearSpringForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "meineke_spring_system.arch";

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");

            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

            std::vector<TissueCell> cells;
            TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                cell.SetBirthTime(-50.0);
                cells.push_back(cell);
            }

            MeshBasedTissue<2> tissue(mesh,cells);
            LinearSpringWithVariableSpringConstantsForce<2> linear_force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            LinearSpringWithVariableSpringConstantsForce<2>* const p_linear_force = &linear_force;

            p_linear_force->UseCutoffPoint(1.1);
            p_linear_force->SetEdgeBasedSpringConstant(true);
            p_linear_force->SetMutantSprings(true, 0.2, 0.3);
            p_linear_force->SetBetaCateninSprings(true);
            p_linear_force->SetApoptoticSprings(true);

            output_arch << p_linear_force;
        }

        {
            ArchiveLocationInfo::SetMeshPathname("mesh/test/data","square_2_elements");

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            LinearSpringWithVariableSpringConstantsForce<2>* p_linear_force;

            // Restore from the archive
            input_arch >> p_linear_force;

            // Test the member data
            TS_ASSERT_EQUALS(p_linear_force->mUseCutoffPoint, true);
            TS_ASSERT_EQUALS(TissueConfig::Instance()->GetMechanicsCutOffLength(), 1.1);
            TS_ASSERT_EQUALS(p_linear_force->mUseEdgeBasedSpringConstant, true);
            TS_ASSERT_EQUALS(p_linear_force->mUseEdgeBasedSpringConstant, true);
            TS_ASSERT_EQUALS(p_linear_force->mUseMutantSprings, true);
            TS_ASSERT_EQUALS(p_linear_force->mUseBCatSprings, true);
            TS_ASSERT_EQUALS(p_linear_force->mUseApoptoticSprings, true);
            TS_ASSERT_DELTA(p_linear_force->mMutantMutantMultiplier, 0.2, 1e-12);
            TS_ASSERT_DELTA(p_linear_force->mNormalMutantMultiplier, 0.3, 1e-12);

            delete p_linear_force;
        }
    }
};

#endif /*TESTFORCES_HPP_*/

