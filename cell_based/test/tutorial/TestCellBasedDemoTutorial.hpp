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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTCELLBASEDDEMOTUTORIAL_HPP_
#define TESTCELLBASEDDEMOTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and cell-based simulations in Chaste =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualise various
 * cell-based simulations. We begin with a simple monolayer simulations, see how to
 *   * change the cell-level model;
 *   * how to impose boundaries;
 *   * how to impose periodic conditions; and
 *   * how to change cell-cycle models.
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * We begin by including the necessary header files. These will be described in detail in 
 * subsequent cell-based tutorials.
 */
#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "RepulsionForce.hpp"
#include "NagaiHondaForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "SmartPointers.hpp"

/*
 * Next, we define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}
 * and defines some test methods.
 */
class TestCellBasedDemoTutorial : public AbstractCellBasedTestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1 - a basic vertex-based simulation ==
     *
     * EMPTYLINE
     *
     * In the first test, we run a simple vertex-based simulation of an epithelial monolayer. 
     * Each cell in the simulation is assigned a simple stochastic cell-cycle model.
     */
    void TestVertexBasedMonolayer() throw (Exception)
    {
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), TRANSIT);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo1");
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(2.0);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/CellBasedDemo1/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
     * java executable.
     *
     * EMPTYLINE
     *
     * == Test 2 - basic node-based simulation ==
     *
     * EMPTYLINE
     *
     * We next show how to modify the previous test to implement a 'node-based' simulation, 
     * in which cells are represented by overlapping spheres (actually circles, since we're 
     * in 2D).
     */
    void TestNodeBasedMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2); //**Changed**//
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //**Changed**//
        NodesOnlyMesh<2> mesh; //**Changed**//
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh); //**Changed**//

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), TRANSIT);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);//**Changed**//
        cell_population.SetMechanicsCutOffLength(1.5); //**Changed**// // to speed up simulations

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo2"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(4.0); //**Changed**//

        MAKE_PTR(RepulsionForce<2>, p_force); //**Changed**//
        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(RandomCellKiller<2>, p_cell_killer, (&cell_population, 0.01)); //**Changed**//
        simulator.AddCellKiller(p_cell_killer);

        simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/CellBasedDemo2/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * EMPTYLINE
     *
     * == Test 3 - basic mesh-based simulation ==
     *
     * We next show how to modify the previous test to implement a 'mesh-based' simulation, 
     * in which cells are represented by their centres and a Voronoi tessellation is used to 
     * find nearest neighbours.
     */
    void TestMeshBasedMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), TRANSIT);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells); //**Changed**//

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo3"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(4.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force); //**Changed**//
        simulator.AddForce(p_force);

        simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     * The results may be visualized using {{{Visualize2dCentreCells}}} as described in the 
     * previous test, with the results directory changed from {{{CellBasedDemo2}}} to {{{CellBasedDemo3}}}.
     *
     * EMPTYLINE
     *
     * == Test 4 - basic mesh-based simulation with ghost nodes ==
     *
     * We next show how to modify the previous test to include 'ghost nodes', which do not 
     * correspond to cells but are sometimes needed when using a Voronoi tessellation. We 
     * will discuss ghost nodes in more detail in subsequent cell-based tutorials.
     */
    void TestMeshBasedMonolayerWithGhostNodes() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2, 2); //**Changed**//
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();//**Changed**//
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT); //**Changed**//

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); //**Changed**//

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo4"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(4.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     * The results may be visualized using {{{Visualize2dCentreCells}}} as described in the 
     * previous test, with the results directory changed from {{{CellBasedDemo3}}} to {{{CellBasedDemo4}}}.
     *
     * EMPTYLINE
     *
     * == Test 5 - basic periodic mesh-based simulation ==
     * 
     * We next show how to modify the previous test to implement a periodic boundary to the 
     * left and right of the domain.
     */
    void TestMeshBasedMonolayerPeriodic() throw (Exception)
    {
        CylindricalHoneycombMeshGenerator generator(5, 2, 2); //**Changed**//
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh(); //**Changed**//

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo5"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(4.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     * The results may be visualized using {{{Visualize2dCentreCells}}} as described in the 
     * previous test, with the results directory changed from {{{CellBasedDemo4}}} to {{{CellBasedDemo5}}}.
     *
     * EMPTYLINE
     *
     * == Test 6 - basic periodic mesh-based simulation with obstructions ==
     *
     * In the final test, we show how to modify the previous test to include one 
     * or more 'obstructions' within the domain. '''THIS IS NOT YET IMPLEMENTED'''
     */
    void TestMeshBasedMonolayerPeriodicSolidBottomBoundary() throw (Exception)
    {
        CylindricalHoneycombMeshGenerator generator(5, 2, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        OffLatticeSimulation<2> simulator(cell_population); //**Changed**//
        simulator.SetOutputDirectory("CellBasedDemo6"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(4.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();
    }
    /*
     * EMPTYLINE
     * 
     * The results may be visualized using {{{Visualize2dCentreCells}}} as described in the 
     * previous test, with the results directory changed from {{{CellBasedDemo5}}} to {{{CellBasedDemo6}}}.
     *
     * EMPTYLINE
     */
};

#endif /*TESTCELLBASEDDEMOTUTORIAL_HPP_*/
