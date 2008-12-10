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
#ifndef TESTFORCES_HPP_
#define TESTFORCES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "FixedCellCycleModelCellsGenerator.hpp"
#include "IngeWntSwatCellCycleModelCellsGenerator.hpp"
#include "MeshBasedTissueWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeinekeInteractionForce.hpp"
#include "ChemotacticForce.hpp"
#include "CellwiseDataGradient.hpp"
#include "WntConcentration.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestForces : public AbstractCancerTestSuite
{
public:

    void TestMeinekeInteractionForceMethods() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();

        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellMutationState mutation_state = HEALTHY;
            if (i==60)
            {
                mutation_state = APC_TWO_HIT;
            }

            TissueCell cell(STEM, mutation_state, new FixedCellCycleModel());
            cell.SetLocationIndex(i);
            cell.SetBirthTime(-10);
            cells.push_back(cell);
        }

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);        
        MeinekeInteractionForce<2> meineke_force;

        /*
         ************************************************************************
         ************************************************************************
         *  Test node velocity calculation
         ************************************************************************
         ************************************************************************
         */

        // Initialise a vector of node velocities
        std::vector<c_vector<double, 2> > node_velocities;
        node_velocities.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_velocities.push_back(zero_vector<double>(2));
        }
        
        // Add velocity contribution from MeinekeInteractionForce 
        /// \todo eventually this should be a force contribution (see #627)
        meineke_force.AddVelocityContribution(node_velocities, tissue);

        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            std::set<unsigned>::iterator iter = ghost_node_indices.find(i);
            bool is_a_ghost_node = (iter!=ghost_node_indices.end());

            if (!is_a_ghost_node)
            {
                TS_ASSERT_DELTA(node_velocities[i][0], 0.0, 1e-4);
                TS_ASSERT_DELTA(node_velocities[i][1], 0.0, 1e-4);
            }
        }

        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point = p_mesh->GetNode(59)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];

        p_mesh->SetNode(59, new_point, false);
        
        // Initialise a vector of new node velocities
        std::vector<c_vector<double, 2> > new_node_velocities;
        new_node_velocities.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             new_node_velocities.push_back(zero_vector<double>(2));
        }
        meineke_force.AddVelocityContribution(new_node_velocities, tissue);

        TS_ASSERT_DELTA(new_node_velocities[60][0], 0.5*p_params->GetSpringStiffness()/p_params->GetDampingConstantMutant(), 1e-4);
        TS_ASSERT_DELTA(new_node_velocities[60][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_velocities[59][0], (-3+4.0/sqrt(7))*p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal(), 1e-4);
        TS_ASSERT_DELTA(new_node_velocities[59][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_velocities[58][0], 0.5*p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal(), 1e-4);
        TS_ASSERT_DELTA(new_node_velocities[58][1], 0.0, 1e-4);

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

        force_on_spring = meineke_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
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

        meineke_force.UseCutoffPoint(dist-0.1);

        force_on_spring = meineke_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                   p_element->GetNodeGlobalIndex(0),
                                                                   tissue);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);
    }


    void TestMeinekeInteractionForceWithEdgeLengthBasedSpring() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false, crypt_width/cells_across);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true); // true = mature cells

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);
        MeinekeInteractionForce<2> meineke_force;

        // Check that the force between nodes is correctly calculated when the 'spring constant' is constant
        meineke_force.SetEdgeBasedSpringConstant(false);

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
            spring_iterator != tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = meineke_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                 nodeB_global_index,
                                                                                 tissue);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],6.25,1e-3);
        }

        // Check that the force between nodes is correctly calculated when the 'spring constant'
        // is proportional to the length of the edge between adjacent cells
        meineke_force.SetEdgeBasedSpringConstant(true);
        tissue.CreateVoronoiTessellation();  // normally done in a simulation loop

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
            spring_iterator != tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = meineke_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                 nodeB_global_index,
                                                                                 tissue);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],4.34027778,1e-3);
        }

        // Choose two interior neighbour nodes
        c_vector<double, 2> force = meineke_force.CalculateForceBetweenNodes(20u, 21u, tissue);
        TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],4.34027778,1e-3);

        // Now move node 21 a bit and check that the force calculation changes correctly
        c_vector<double,2> shift;
        shift[0] = 0.01;
        shift[1] = 0.0;
        ChastePoint<2> new_point(p_mesh->GetNode(21u)->rGetLocation() + shift);
        p_mesh->SetNode(21u, new_point, false);

        // Check that the new force between nodes is correctly calculated
        tissue.CreateVoronoiTessellation();
        c_vector<double, 2> new_force = meineke_force.CalculateForceBetweenNodes(20u, 21u, tissue);

        // Force calculation: shift is along x-axis so we should have
        // new_edge_length = (5/6 + shift[0])*tan(0.5*arctan(5*sqrt(3)/(5 + 12*shift[0]))),
        // force^2 = mu^2 * (new_edge_length*sqrt(3))^2 * (1 - 5/6 - shift[0])^2
        TS_ASSERT_DELTA(new_force[0]*new_force[0] + new_force[1]*new_force[1], 3.83479824,1e-3);
    }

    // Test on a periodic mesh
    void TestMeinekeInteractionForceWithEdgeBasedSpringsOnPeriodicMesh() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);// true = mature cells

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);
        MeinekeInteractionForce<2> meineke_force;

        // Check that the force between nodes is correctly calculated when the spring constant is constant (!)
        meineke_force.SetEdgeBasedSpringConstant(false);

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
            spring_iterator != tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = meineke_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                 nodeB_global_index,
                                                                                 tissue);

            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 6.25, 1e-3);
        }

        // Check that the force between nodes is correctly calculated when the spring constant
        // is proportional to the length of the edge between adjacenet cells
        meineke_force.SetEdgeBasedSpringConstant(true);
        tissue.CreateVoronoiTessellation();

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
            spring_iterator != tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = meineke_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                 nodeB_global_index,
                                                                                 tissue);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 4.34027778, 1e-3);
        }
    }


    void TestMeinekeInteractionForceWithAreaBasedVisocity() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 2;

        // Test a non-periodic mesh
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Scale the mesh in one direction
        p_mesh->Scale(1.0,1.2);

        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);// true = mature cells

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        crypt.CreateVoronoiTessellation();  // normally done in a simulation loop

        MeinekeInteractionForce<2> meineke_force;
        
        // Initialise a vector of node velocities
        std::vector<c_vector<double, 2> > node_velocities;
        node_velocities.reserve(crypt.GetNumNodes());
        
        for (unsigned i=0; i<crypt.GetNumNodes(); i++)
        {
             node_velocities.push_back(zero_vector<double>(2));
        }
        
        // Add velocity contribution from MeinekeInteractionForce 
        /// \todo eventually this should be a force contribution (see #627)
        meineke_force.AddVelocityContribution(node_velocities, crypt);
        
        std::vector<double> norm_vel;

        for (unsigned i=0; i<node_velocities.size(); i++)
        {
            // Check if this is a real cell
            if (ghost_node_indices.find(i)==ghost_node_indices.end())
            {
                norm_vel.push_back(norm_2(node_velocities[i]));
            }
        }

        // Now check that the velocities scale correctly when the viscosity is area-dependent
        meineke_force.SetAreaBasedViscosity(true);
        
        // Initialise a new vector of node velocities
        std::vector<c_vector<double, 2> > new_node_velocities;
        new_node_velocities.reserve(crypt.GetNumNodes());
        
        for (unsigned i=0; i<crypt.GetNumNodes(); i++)
        {
             new_node_velocities.push_back(zero_vector<double>(2));
        }
        
        // Add velocity contribution from MeinekeInteractionForce 
        /// \todo eventually this should be a force contribution (see #627)
        meineke_force.AddVelocityContribution(new_node_velocities, crypt);
        
        std::vector<double> norm_vel_area;

        for (unsigned i=0; i<new_node_velocities.size(); i++)
        {
            // Check if this is a real cell
            if (ghost_node_indices.find(i)==ghost_node_indices.end())
            {
                norm_vel_area.push_back(norm_2(new_node_velocities[i]));
            }
        }

        TS_ASSERT(norm_vel.size() > 0);

        // Note that d0 and d1 are hardcoded in TissueSimulation::mpMechanicsSystem->rCalculateVelocitiesOfEachNode()
        for (unsigned i=0; i<norm_vel.size(); i++)
        {
            TS_ASSERT_DELTA(norm_vel_area[i], norm_vel[i]/(0.1 +  1.2*0.9), 1e-3);
        }
    }


    void TestMeinekeInteractionForceWithAreaBasedVisocityOnAPeriodicMesh() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 2;

        double scale_factor = 1.2;
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, scale_factor);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true); // true = mature cells

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);

        MeinekeInteractionForce<2> meineke_force;

        // It seems quite difficult to test this on a periodic mesh,
        // so just check the areas of all the cells are correct.
        
        /// \todo this doesn't test any method on MeinekeInteractionForce, 
        // so perhaps should moved to TestMeshBasedTissue?
        tissue.CreateVoronoiTessellation();
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            // Check if this is a real cell
            if (ghost_node_indices.find(i)==ghost_node_indices.end())
            {
                double area = tissue.rGetVoronoiTessellation().GetFaceArea(i);
                TS_ASSERT_DELTA(area, sqrt(3)*scale_factor*scale_factor/2, 1e-6);
            }
        }
    }


    void TestMeinekeInteractionForceWithSpringConstantsForMutantCells()
    {
        // Create a small tissue
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0));
        nodes.push_back(new Node<2>(1, false, 0, 2));
        nodes.push_back(new Node<2>(2, false, 2, 2));
        nodes.push_back(new Node<2>(3, false, 2, 0));

        MutableMesh<2,2> mesh(nodes);

        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh);

        MeshBasedTissue<2> tissue(mesh, cells);

        MeinekeInteractionForce<2> meineke_force;

        // Set cells mutation states
        tissue.rGetCellUsingLocationIndex(0).SetMutationState(HEALTHY);
        tissue.rGetCellUsingLocationIndex(1).SetMutationState(LABELLED);
        tissue.rGetCellUsingLocationIndex(2).SetMutationState(APC_TWO_HIT);
        tissue.rGetCellUsingLocationIndex(3).SetMutationState(BETA_CATENIN_ONE_HIT);

        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(0, 1, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(1, 2, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(2, 3, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(3, 0, tissue)), 15.0, 1e-10);

        meineke_force.SetMutantSprings(true);

        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(0, 1, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(1, 2, tissue)), 22.5, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(2, 3, tissue)), 30.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(3, 0, tissue)), 22.5, 1e-10);

        meineke_force.SetMutantSprings(true, 4.0, 3.0);

        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(0, 1, tissue)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(1, 2, tissue)), 45.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(2, 3, tissue)), 60.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(3, 0, tissue)), 45.0, 1e-10);
    }


    void TestMeinekeInteractionForceWithSpringConstantsForIngeBCatCells()
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(6, 12, 0, true, 1.1);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        IngeWntSwatCellCycleModelCellsGenerator<2> cells_generator(2u);
        cells_generator.GenerateForCrypt(cells, *p_mesh, false);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        // As there is no tissue simulation, we must explicitly initialise the cells
        crypt.InitialiseCells();

        MeinekeInteractionForce<2> meineke_force;

        TS_ASSERT_DELTA(norm_2(meineke_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.50, 1e-10);

        meineke_force.SetBCatSprings(true);
        crypt.CreateVoronoiTessellation();  // this method is normally called in a simulation loop

        /// \todo this is currently a rather poor test - it just checks that
        /// there is SOME dependency of the spring constant on the betat catenin level
        /// experienced by both cells
        
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.5*8.59312/18.14, 1e-5);
        
        CancerParameters::Instance()->SetBetaCatSpringScaler(20/6.0);        
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.5*8.59312/20.0, 1e-5);

        WntConcentration::Destroy();
    }

    void TestMeinekeInteractionForceWithSpringConstantsForApoptoticCells()
    {
        // Set up stretched tissue
        HoneycombMeshGenerator generator(4, 4, 0, false, 2.0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, *p_mesh);

        MeshBasedTissueWithGhostNodes<2> stretched_tissue(*p_mesh, cells, ghost_node_indices);

        // As there is no tissue simulation we must explicitly initialise the cells
        stretched_tissue.InitialiseCells();

        // Set one of the non-boundary cells to be necrotic
        stretched_tissue.rGetCellUsingLocationIndex(6).SetCellType(APOPTOTIC);

        MeinekeInteractionForce<2> meineke_force;
        meineke_force.SetApoptoticSprings(true);

        TS_ASSERT_EQUALS( stretched_tissue.rGetCellUsingLocationIndex(6).GetCellType(), APOPTOTIC);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(6, 10, stretched_tissue)), 3.3333, 1e-4);

        // Set a neighbouring cell to be necrotic
        stretched_tissue.rGetCellUsingLocationIndex(10).SetCellType(APOPTOTIC);

        TS_ASSERT_EQUALS( stretched_tissue.rGetCellUsingLocationIndex(10).GetCellType(), APOPTOTIC);
        TS_ASSERT_DELTA( norm_2(meineke_force.CalculateForceBetweenNodes(6, 10, stretched_tissue)), 1.8750, 1e-4);

        // Now do similar tests for a squashed tissue
        HoneycombMeshGenerator generator2(4, 4, 0, false, 0.5);
        MutableMesh<2,2>* p_mesh2 = generator2.GetMesh();
        std::set<unsigned> ghost_node_indices2 = generator2.GetGhostNodeIndices();

        std::vector<TissueCell> cells2;
        FixedCellCycleModelCellsGenerator<2> cells_generator2;
        cells_generator2.GenerateBasic(cells2, *p_mesh2);

        MeshBasedTissueWithGhostNodes<2> squashed_tissue(*p_mesh2, cells2, ghost_node_indices2);
        squashed_tissue.InitialiseCells();

        squashed_tissue.rGetCellUsingLocationIndex(6).SetCellType(APOPTOTIC);

        MeinekeInteractionForce<2> meineke_force2;
        meineke_force2.SetApoptoticSprings(true);

        TS_ASSERT_DELTA( norm_2(meineke_force2.CalculateForceBetweenNodes(6, 10, squashed_tissue)), 4.0909, 1e-4);

        squashed_tissue.rGetCellUsingLocationIndex(10).SetCellType(APOPTOTIC);

        TS_ASSERT_DELTA( norm_2(meineke_force2.CalculateForceBetweenNodes(6, 10, squashed_tissue)), 2.8125, 1e-4);
    }


    void TestMeinekeInteractionForceCalculationIn3d() throw (Exception)
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;
        TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            cell.SetLocationIndex(i);
            cell.SetBirthTime(-50.0);
            cells.push_back(cell);
        }

        MeshBasedTissue<3> tissue(mesh,cells);
        MeinekeInteractionForce<3> meineke_force;

        // Test forces on springs
        unsigned nodeA = 0, nodeB = 1;
        Element<3,3>* p_element = mesh.GetElement(0);
        c_vector<double, 3> force = meineke_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(nodeA),
                                                                             p_element->GetNodeGlobalIndex(nodeB),
                                                                             tissue);
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(force[i], 0.0, 1e-6);
        }

        // Initialise a vector of node velocities
        std::vector<c_vector<double, 3> > node_velocities;
        node_velocities.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_velocities.push_back(zero_vector<double>(3));
        }
        
        // Add velocity contribution from MeinekeInteractionForce 
        /// \todo eventually this should be a force contribution (see #627)
        meineke_force.AddVelocityContribution(node_velocities, tissue);

        for (unsigned j=0; j<4; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(node_velocities[j](k), 0.0, 1e-6);
            }
        }

        // Scale entire mesh and check that forces are correctly calculated
        CancerParameters *p_params = CancerParameters::Instance();
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
        
        // Recalculate node velocities (we can just re-use node_velocities,
        // as previously each node had zero velocity)
        meineke_force.AddVelocityContribution(node_velocities, tissue);

        for (unsigned j=0; j<4; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(fabs(node_velocities[j](k)), p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(scale_factor-1)*sqrt(2),1e-6);
            }
        }

        // Move one node and check that forces are correctly calculated
        MutableMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        
        MeshBasedTissue<3> tissue2(mesh2,cells);
        MeinekeInteractionForce<3> meineke_force2;
        
        c_vector<double,3> old_point = mesh2.GetNode(0)->rGetLocation();
        ChastePoint<3> new_point;
        new_point.rGetLocation()[0] = 0.0;
        new_point.rGetLocation()[1] = 0.0;
        new_point.rGetLocation()[2] = 0.0;
        mesh2.SetNode(0, new_point, false);

        unsigned nodeA2 = 0, nodeB2 = 1;
        Element<3,3>* p_element2 = mesh2.GetElement(0);
        c_vector<double,3> force2 = meineke_force2.CalculateForceBetweenNodes(p_element2->GetNodeGlobalIndex(nodeA2),
                                                                              p_element2->GetNodeGlobalIndex(nodeB2),
                                                                              tissue2);

        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(fabs(force2[i]),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }

        // Initialise a vector of node velocities
        std::vector<c_vector<double,3> > node_velocities2;
        node_velocities2.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_velocities2.push_back(zero_vector<double>(3));
        }
        
        // Add velocity contribution from MeinekeInteractionForce 
        /// \todo eventually this should be a force contribution (see #627)
        meineke_force2.AddVelocityContribution(node_velocities2, tissue2);

        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(node_velocities2[0](i),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
    }


    void TestMeinekeInteractionForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "meineke_spring_system.arch";

        unsigned num_nodes;
        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");

            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            num_nodes = mesh.GetNumNodes();

            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

            std::vector<TissueCell> cells;
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                cell.SetLocationIndex(i);
                cell.SetBirthTime(-50.0);
                cells.push_back(cell);
            }

            MeshBasedTissue<2> tissue(mesh,cells);
            MeinekeInteractionForce<2> meineke_force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            MeinekeInteractionForce<2>* const p_meineke_force = &meineke_force;
            
            p_meineke_force->UseCutoffPoint(1.1);
            p_meineke_force->SetAreaBasedViscosity(true);
            p_meineke_force->SetEdgeBasedSpringConstant(true);
            p_meineke_force->SetMutantSprings(true,0.2,0.3);
            p_meineke_force->SetBCatSprings(true);
            p_meineke_force->SetApoptoticSprings(true);

            output_arch << p_meineke_force;
        }

        {
            MeshBasedTissue<2>::meshPathname = "mesh/test/data/square_2_elements";
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            MeinekeInteractionForce<2>* p_meineke_force;

            // Restore from the archive
            input_arch >> p_meineke_force;

            // Test the member data
            TS_ASSERT_EQUALS(p_meineke_force->mUseCutoffPoint,true);            
            TS_ASSERT_EQUALS(p_meineke_force->mUseEdgeBasedSpringConstant, true);
            TS_ASSERT_EQUALS(p_meineke_force->mUseAreaBasedViscosity, true);
            TS_ASSERT_EQUALS(p_meineke_force->mUseEdgeBasedSpringConstant, true);
            TS_ASSERT_EQUALS(p_meineke_force->mUseMutantSprings, true);
            TS_ASSERT_EQUALS(p_meineke_force->mUseBCatSprings, true);
            TS_ASSERT_EQUALS(p_meineke_force->mUseApoptoticSprings, true);
            TS_ASSERT_DELTA(p_meineke_force->mCutoffPoint,1.1,1e-12);
            TS_ASSERT_DELTA(p_meineke_force->mMutantMutantMultiplier, 0.2, 1e-12);
            TS_ASSERT_DELTA(p_meineke_force->mNormalMutantMultiplier, 0.3, 1e-12);

            delete p_meineke_force;
        }
    }
    
    void TestChemotacticForceMethods() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 0;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, LABELLED, new FixedCellCycleModel());
            cell.SetLocationIndex(i);
            cell.SetBirthTime(-10);
            cells.push_back(cell);
        }

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);

        // Set up cellwise data and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(), 1);
        p_data->SetTissue(tissue);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            p_data->SetValue(x/50.0, p_mesh->GetNode(i));
        }

        ChemotacticForce<2> chemotactic_force;
        
        // Initialise a vector of new node velocities
        std::vector<c_vector<double, 2> > node_velocities;
        node_velocities.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_velocities.push_back(zero_vector<double>(2));
        }
        chemotactic_force.AddVelocityContribution(node_velocities, tissue);

        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            bool is_a_ghost_node = tissue.rGetGhostNodes()[i];

            if (!is_a_ghost_node)
            {
                double x = p_mesh->GetNode(i)->rGetLocation()[0];
                double c = x/50;
                double norm_grad_c = 1.0/50.0;
                double force_magnitude = chemotactic_force.GetChemotacticForceMagnitude(c, norm_grad_c);

                // As only labelled cells experience the chemotactic force, we must use
                // the mutant damping constant
                double damping = CancerParameters::Instance()->GetDampingConstantMutant();

                // Fc = force_magnitude*(1,0), Fspring=0 => velocity = damping*force_magnitude*(1,0)
                TS_ASSERT_DELTA(node_velocities[i][0], force_magnitude/damping, 1e-4);
                TS_ASSERT_DELTA(node_velocities[i][1], 0.0, 1e-4);
            }
        }
    }

    void TestChemotacticForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "chemotaxis_spring_system.arch";

        unsigned num_nodes;
        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");

            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            num_nodes = mesh.GetNumNodes();

            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

            std::vector<TissueCell> cells;
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                cell.SetLocationIndex(i);
                cell.SetBirthTime(-50.0);
                cells.push_back(cell);
            }

            MeshBasedTissue<2> tissue(mesh,cells);
            ChemotacticForce<2> chemotactic_force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            ChemotacticForce<2> * const p_chemotactic_force = &chemotactic_force;

            p_chemotactic_force->SetAreaBasedViscosity(true);

            output_arch << p_chemotactic_force;
        }

        {
            MeshBasedTissue<2>::meshPathname = "mesh/test/data/square_2_elements";
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            ChemotacticForce<2>* p_chemotactic_force;

            // Restore from the archive
            input_arch >> p_chemotactic_force;

            // Test the member data
            TS_ASSERT_EQUALS(p_chemotactic_force->mUseAreaBasedViscosity, true);

            delete p_chemotactic_force;
        }
    }
    
};

#endif /*TESTFORCES_HPP_*/

