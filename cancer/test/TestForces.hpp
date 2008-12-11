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
#include "MeinekeInteractionWithVariableSpringConstantsForce.hpp"
#include "ChemotacticForce.hpp"
#include "CellwiseDataGradient.hpp"
#include "CryptProjectionForce.hpp"
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

        meineke_force.AddForceContribution(node_forces, tissue);

        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            std::set<unsigned>::iterator iter = ghost_node_indices.find(i);
            bool is_a_ghost_node = (iter!=ghost_node_indices.end());

            if (!is_a_ghost_node)
            {
                TS_ASSERT_DELTA(node_forces[i][0], 0.0, 1e-4);
                TS_ASSERT_DELTA(node_forces[i][1], 0.0, 1e-4);
            }
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
        meineke_force.AddForceContribution(new_node_forces, tissue);

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
        MeinekeInteractionWithVariableSpringConstantsForce<2> meineke_force;

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
        MeinekeInteractionWithVariableSpringConstantsForce<2> meineke_force;

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

        MeinekeInteractionWithVariableSpringConstantsForce<2> meineke_force;

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

        MeinekeInteractionWithVariableSpringConstantsForce<2> meineke_force;

        TS_ASSERT_DELTA(norm_2(meineke_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.50, 1e-10);

        meineke_force.SetBCatSprings(true);
        crypt.CreateVoronoiTessellation();  // this method is normally called in a simulation loop

        /// \todo this is currently a rather poor test - it just checks that
        /// there is SOME dependency of the spring constant on the beta catenin level
        /// experienced by both cells (see #627)
        
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

        MeinekeInteractionWithVariableSpringConstantsForce<2> meineke_force;
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

        MeinekeInteractionWithVariableSpringConstantsForce<2> meineke_force2;
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

        // Initialise a vector of node forces
        std::vector<c_vector<double, 3> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(3));
        }

        meineke_force.AddForceContribution(node_forces, tissue);

        for (unsigned j=0; j<4; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(node_forces[j](k), 0.0, 1e-6);
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
        
        // Recalculate node forces (we can just re-use node_forces,
        // as previously each node had zero net force on it)
        meineke_force.AddForceContribution(node_forces, tissue);

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
            TS_ASSERT_DELTA(fabs(force2[i]),p_params->GetSpringStiffness()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }

        // Initialise a vector of node forces
        std::vector<c_vector<double,3> > node_forces2;
        node_forces2.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces2.push_back(zero_vector<double>(3));
        }
        
        meineke_force2.AddForceContribution(node_forces2, tissue2);

        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(node_forces2[0](i),p_params->GetSpringStiffness()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
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
            MeinekeInteractionWithVariableSpringConstantsForce<2> meineke_force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            MeinekeInteractionWithVariableSpringConstantsForce<2>* const p_meineke_force = &meineke_force;
            
            p_meineke_force->UseCutoffPoint(1.1);
            p_meineke_force->SetEdgeBasedSpringConstant(true);
            p_meineke_force->SetMutantSprings(true, 0.2, 0.3);
            p_meineke_force->SetBCatSprings(true);
            p_meineke_force->SetApoptoticSprings(true);

            output_arch << p_meineke_force;
        }

        {
            MeshBasedTissue<2>::meshPathname = "mesh/test/data/square_2_elements";
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            MeinekeInteractionWithVariableSpringConstantsForce<2>* p_meineke_force;

            // Restore from the archive
            input_arch >> p_meineke_force;

            // Test the member data
            TS_ASSERT_EQUALS(p_meineke_force->mUseCutoffPoint,true);            
            TS_ASSERT_EQUALS(p_meineke_force->mUseEdgeBasedSpringConstant, true);
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
        
        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }
        chemotactic_force.AddForceContribution(node_forces, tissue);

        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            bool is_a_ghost_node = tissue.rGetGhostNodes()[i];

            if (!is_a_ghost_node)
            {
                double x = p_mesh->GetNode(i)->rGetLocation()[0];
                double c = x/50;
                double norm_grad_c = 1.0/50.0;
                double force_magnitude = chemotactic_force.GetChemotacticForceMagnitude(c, norm_grad_c);
                
                // Fc = force_magnitude*(1,0), Fspring=0
                TS_ASSERT_DELTA(node_forces[i][0], force_magnitude, 1e-4);
                TS_ASSERT_DELTA(node_forces[i][1], 0.0, 1e-4);
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

            /// \todo Test the member data (see #627)
            delete p_chemotactic_force;
        }
    }
  
    void TestCryptProjectionForceMethods() throw (Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();

        // Create a mesh
        unsigned num_cells_width = 10;
        unsigned num_cells_depth = 10;
        unsigned thickness_of_ghost_layer = 0;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Centre the mesh at (0,0)
        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);

        double width_of_mesh = width_extremes[1] - width_extremes[0];
        double height_of_mesh = height_extremes[1] - height_extremes[0];

        p_mesh->Translate(-width_of_mesh/2, -height_of_mesh/2);

        // Create some cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            cell.SetLocationIndex(i);

            if (i==4 || i==5)
            {
                cell.SetBirthTime(-0.5);
            }
            else
            {
                cell.SetBirthTime(-10.0);
            }
            cells.push_back(cell);
        }

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        tissue.MarkSpring(tissue.rGetCellUsingLocationIndex(4), tissue.rGetCellUsingLocationIndex(5));

        // Create a spring system with crypt surface z = 2*r
        p_params->SetCryptProjectionParameterA(2.0);
        p_params->SetCryptProjectionParameterB(1.0);       
        CryptProjectionForce crypt_projection_force;

        // Test get methods
        TS_ASSERT_DELTA(crypt_projection_force.GetA(), 2.0, 1e-12);
        TS_ASSERT_DELTA(crypt_projection_force.GetB(), 1.0, 1e-12);

        // Test crypt height and gradient calculations
        c_vector<double, 2> node_location_2d = p_mesh->GetNode(0)->rGetLocation();
        TS_ASSERT_DELTA(crypt_projection_force.CalculateCryptSurfaceHeightAtPoint(node_location_2d), 2.0*pow(norm_2(node_location_2d),1.0), 1e-12);
        TS_ASSERT_DELTA(crypt_projection_force.CalculateCryptSurfaceDerivativeAtPoint(node_location_2d), 2.0, 1e-12);

        // Test updating of mNode3dLocationMap
        crypt_projection_force.UpdateNode3dLocationMap(tissue);

        // Move a node slightly
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = node_location_2d[0]+0.05;
        new_point.rGetLocation()[1] = node_location_2d[1];
        p_mesh->SetNode(0, new_point, false);
 

        // Test UpdateNode3dLocationMap()

        c_vector<double, 2> new_node_location_2d;
        new_node_location_2d[0] = new_point.rGetLocation()[0];
        new_node_location_2d[1] = new_point.rGetLocation()[1];

        crypt_projection_force.UpdateNode3dLocationMap(tissue);

        // Check the map updates correctly (note that we have used no ghost nodes, so the map does contain 0)
        c_vector<double, 3> calculated_new_node_location_3d = crypt_projection_force.mNode3dLocationMap[0];
        c_vector<double, 3> correct_new_node_location_3d;

        correct_new_node_location_3d[0] = new_node_location_2d[0];
        correct_new_node_location_3d[1] = new_node_location_2d[1];
        correct_new_node_location_3d[2] = crypt_projection_force.CalculateCryptSurfaceHeightAtPoint(new_node_location_2d);

        TS_ASSERT_DELTA(calculated_new_node_location_3d[0], correct_new_node_location_3d[0], 1e-12);
        TS_ASSERT_DELTA(calculated_new_node_location_3d[1], correct_new_node_location_3d[1], 1e-12);
        TS_ASSERT_DELTA(calculated_new_node_location_3d[2], correct_new_node_location_3d[2], 1e-12);

   
        // Test force calculation on a normal spring

        c_vector<double,2> force_on_spring; // between nodes 0 and 1

        // Find one of the elements that nodes 0 and 1 live on
        ChastePoint<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0] + 0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01;

        unsigned elem_index = p_mesh->GetContainingElementIndex(new_point2, false);
        Element<2,2>* p_element = p_mesh->GetElement(elem_index);

        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                            p_element->GetNodeGlobalIndex(0),
                                                                            tissue);

        TS_ASSERT_DELTA(force_on_spring[0], -5.7594, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0230 , 1e-4);

///////////////////////////////////////////////////////////////// 
        // Test force calculation with a cutoff

        double dist = norm_2( p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(), 
                              p_element->GetNode(1)->rGetLocation()) );
                              
        crypt_projection_force.UseCutoffPoint(dist - 0.1);

        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                            p_element->GetNodeGlobalIndex(0),
                                                                            tissue);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);


        // Test force calculation for a pair of newly born neighbouring cells
        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(4, 5, tissue);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);

        tissue.UnmarkSpring(tissue.rGetCellUsingLocationIndex(4), tissue.rGetCellUsingLocationIndex(5));

        // Test force calculation for a particular node
        
        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }
        
        crypt_projection_force.AddForceContribution(node_forces, tissue);      

        TS_ASSERT_DELTA(node_forces[0][0], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_forces[0][1], 0.0, 1e-4);

        // Test that in the case of a flat crypt surface (mA=mB=0), the results are the same as for Meineke2001SpringSystem
        p_params->SetCryptProjectionParameterA(0.001);
        p_params->SetCryptProjectionParameterB(0.001);
        CryptProjectionForce flat_crypt_projection_force;
        MeinekeInteractionForce<2> meineke_force;

        // Normally this would be set up at the start of rCalculateforcesOfEachNode
        flat_crypt_projection_force.UpdateNode3dLocationMap(tissue);

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
            spring_iterator != tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

            c_vector<double, 2> force_flat = flat_crypt_projection_force.CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, tissue);
            c_vector<double, 2> force_meineke = meineke_force.CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, tissue);

            TS_ASSERT_DELTA( force_flat[0], force_meineke[0], 1e-3);
            TS_ASSERT_DELTA( force_flat[1], force_meineke[1], 1e-3);
        }
    }

    /**
     * Note: WntBasedChemotaxis should be possible in other spring systems. If/when
     * this is implemented, this test should be moved to somewhere more appropriate.
     */
    void TestCryptProjectionForceWithWntBasedChemotaxis() throw (Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();

        // Create a mesh
        unsigned num_cells_width = 10;
        unsigned num_cells_depth = 10;
        unsigned thickness_of_ghost_layer = 0;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Centre the mesh at (0,0)
        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);

        double width_of_mesh = (num_cells_width/num_cells_width)*(width_extremes[1] - width_extremes[0]);
        double height_of_mesh = (num_cells_depth/num_cells_depth)*(height_extremes[1] - height_extremes[0]);

        p_mesh->Translate(-width_of_mesh/2, -height_of_mesh/2);

        // Create some cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            cell.SetLocationIndex(i);
            cell.SetBirthTime(-10.0);
            cells.push_back(cell);
        }

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        tissue.MarkSpring(tissue.rGetCellUsingLocationIndex(4), tissue.rGetCellUsingLocationIndex(5));

        WntConcentration::Instance()->SetType(RADIAL);
        WntConcentration::Instance()->SetTissue(tissue);

        // Create a spring system with crypt surface z = 2*r
        p_params->SetCryptProjectionParameterA(2.0);
        p_params->SetCryptProjectionParameterB(1.0);
        CryptProjectionForce crypt_projection_force;

        crypt_projection_force.SetWntChemotaxis(false);
               
        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > old_node_forces;
        old_node_forces.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             old_node_forces.push_back(zero_vector<double>(2));
        }
        
        // Calculate node forces
        crypt_projection_force.AddForceContribution(old_node_forces, tissue);
        
        // Store the force of a particular node without Wnt-chemotaxis
        c_vector<double,2> old_force = old_node_forces[11];      

        // Now turn on Wnt-chemotaxis
        crypt_projection_force.SetWntChemotaxis(true);
        
        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > new_node_forces;
        new_node_forces.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             new_node_forces.push_back(zero_vector<double>(2));
        }
                
        // Calculate node forces
        crypt_projection_force.AddForceContribution(new_node_forces, tissue);

        // Store the force of the same node, but now with Wnt-chemotaxis
        c_vector<double,2> new_force = new_node_forces[11];

        double wnt_chemotaxis_strength = CancerParameters::Instance()->GetWntChemotaxisStrength();
        c_vector<double,2> wnt_component = wnt_chemotaxis_strength*WntConcentration::Instance()->GetWntGradient(&(cells[11]));

        TS_ASSERT_DELTA(new_force[0], old_force[0]+wnt_component[0], 1e-4);
        TS_ASSERT_DELTA(new_force[1], old_force[1]+wnt_component[1], 1e-4);
    }

    void TestCryptProjectionForceWithArchiving() throw (Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();

        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "crypt_projection_spring_system.arch";

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

            MeshBasedTissue<2> crypt(mesh,cells);
            p_params->SetCryptProjectionParameterA(1.0);
            p_params->SetCryptProjectionParameterB(2.0);
            CryptProjectionForce crypt_projection_force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            CryptProjectionForce * const p_crypt_projection_force = &crypt_projection_force;

            p_crypt_projection_force->UseCutoffPoint(1.1);

            output_arch << p_crypt_projection_force;
        }

        {
            MeshBasedTissue<2>::meshPathname = "mesh/test/data/square_2_elements";

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            CryptProjectionForce* p_crypt_projection_force;

            // Restore from the archive
            input_arch >> p_crypt_projection_force;

            // Test the member data
            TS_ASSERT_EQUALS(p_crypt_projection_force->mUseCutoffPoint, true);
            TS_ASSERT_DELTA(p_crypt_projection_force->mCutoffPoint, 1.1, 1e-12);
            TS_ASSERT_DELTA(p_crypt_projection_force->GetA(), 1.0, 1e-12);
            TS_ASSERT_DELTA(p_crypt_projection_force->GetB(), 2.0, 1e-12);

            delete p_crypt_projection_force;
        }
    }
  
    void TestForceCollection() throw (Exception)
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

        // Create two different force laws and add to a std::list
        MeinekeInteractionForce<2> meineke_force;
        
        p_params->SetCryptProjectionParameterA(0.0001);
        p_params->SetCryptProjectionParameterB(0.0001);
        CryptProjectionForce crypt_projection_force;
        
        std::vector<AbstractForce<2>* > forces;
        forces.push_back(&meineke_force);
        forces.push_back(&crypt_projection_force);
        
        // Test node force calculation

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());
        
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }
        
        // Add force contributions
        for (std::vector<AbstractForce<2>* >::iterator iter = forces.begin();
             iter != forces.end();
             iter++)
        {
             (*iter)->AddForceContribution(node_forces, tissue);
        }

        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            std::set<unsigned>::iterator iter = ghost_node_indices.find(i);
            bool is_a_ghost_node = (iter!=ghost_node_indices.end());

            if (!is_a_ghost_node)
            {
                TS_ASSERT_DELTA(node_forces[i][0], 0.0, 1e-4);
                TS_ASSERT_DELTA(node_forces[i][1], 0.0, 1e-4);
            }
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

        // Add force contributions
        for (std::vector<AbstractForce<2>* >::iterator iter = forces.begin();
             iter != forces.end();
             iter++)
        {
             (*iter)->AddForceContribution(new_node_forces, tissue);
        }

        // Forces should be twice the forces found using Meineke alone (since a flat crypt is used)
        TS_ASSERT_DELTA(new_node_forces[60][0], 2*0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[60][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[59][0], 2*(-3+4.0/sqrt(7))*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[59][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[58][0], 2*0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[58][1], 0.0, 1e-4);
    }
    
};

#endif /*TESTFORCES_HPP_*/

