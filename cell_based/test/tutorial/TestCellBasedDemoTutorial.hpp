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
 * This tutroial is designed to give you a quick introduction to running cell-based
 * simulations in Chaste. Full details are postponed until later tutorials.
 *
 * We begin with a simple monolayer simulation and see how to:
 *   * change the cell-level model;
 *   * how to impose boundaries;
 *   * how to impose periodic conditions;
 *   * how to specify how to remove cells; and
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
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "AdhesionPottsUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NagaiHondaForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "RepulsionForce.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
/*
 * Next, we define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}
 * and defines some test methods. We are using {{{AbstractCellBasedTestSuite}}} instead of {{{CxxTest::TestSuite}}} as this
 * handles some of the simulations set up for us, details are given in later tutorials,
 * UserTutorials/RunningMeshBasedSimulations and UserTutorials/RunningNodeBasedSimulations.
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
     * Each cell in the simulation is assigned a simple stochastic cell-cycle model, the cells will divide randomly and never stop proliferating.
     */
    void TestVertexBasedMonolayer() throw (Exception)
    {
        /* The first thing we define is a 2D (specified by the <2,2>) mesh which holds the spatial information of the simulation. To do this we use one of a
         * number of {{{MeshGenerators}}}.*/
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* We now generate a collection of cells. We do this by using a {{{CellsGenerator}}} and we specify the proliferative
         * behaviour of the cell by choosing a {{{CellCycleModel}}}, here we choose a {{{StochasticDurationCellCycleModel}}} where
         * each cell is given a division time, drawn from a uniform distribution, when it is created. For a vertex simulation
         * we need as may cells as elements in the mesh.*/
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), TRANSIT);

        /* We now create a {{{CellPopulation}}} object (passing in the mesh and cells) to connect the mesh and the cells together.
         * Here that is a {{{VertexBasedCellPopulation}}} and the dimension is <2>.*/
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We now create an {{{OffSimulation}}} object and pass in the {{{CellPopulation}}}. We also set some
         * options on the simulation like output directory, output multiple (so we don't visualise every timestep),
         * and end time.
         */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo1");
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(20.0);

        /* In order to specify how cells move around we create a "shared pointer" to a
         * {{{Force}}} object and pass it to the {{{OffLatticeSimulation}}}. This is done using the MAKE_PTR macro as follows.
         */
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* Finally we call the {{{Solve}}} method on the simulation to run the simulation.*/
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
        /* We now need to create a {{{NodesOnlyMesh}}} we do this by first creating a {{{MutableMesh}}}
         * and passing this to a helper method {{{ConstructNodesWithoutMesh}}}. Note that we need to
         * delete {{{p_mesh}}} at the end of the test, to avoid memory leaks, as we created this pointer.
         */
        HoneycombMeshGenerator generator(2, 2); //**Changed**//
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //**Changed**//
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>; //**Changed**//
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh); //**Changed**//

        /* We create the cells as before, only this time we need one cell per node.*/
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), TRANSIT); //**Changed**//

        /* This time we create a {{{NodeBasedCellPopulation}}} as we are using a {{{NodesOnlyMesh}}}.
         * We also set a cut off length to speed up simulations. This defines a radius of interaction for the cells.*/
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);//**Changed**//
        cell_population.SetMechanicsCutOffLength(1.5); //**Changed**//

        /* We create an {{{OffSimulation}}} object as before, all we change is the output directory
         * and output results more often as a larger default timestep is used for these simulations. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo2"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(12); //**Changed**//
        simulator.SetEndTime(20.0);

        /* We use a different {{{Force}}} which is suitable for node based simulations.
         */
        MAKE_PTR(RepulsionForce<2>, p_force); //**Changed**//
        simulator.AddForce(p_force);

        /* In all types of simulation you may specify how cells are removed from the simulation by specifying
         * a {{{CellKiller}}}. You create these in the same was as the {{{Force}}} and pass them to the {{{CellBasedSimulation}}}.
         * Note that here the constructor for {{{RandomCellKiller}}} requires some arguments to be passed to it, therefore we use the
         * {{{MAKE_PTR_ARGS}}} macro.
         */
        MAKE_PTR_ARGS(RandomCellKiller<2>, p_cell_killer, (&cell_population, 0.01)); //**Changed**//
        simulator.AddCellKiller(p_cell_killer);

        /* Again we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();

        /* Before finishing the test we delete the pointer to the mesh to avoid memory leaks */
        delete p_mesh; //** Changed**// to stop memory leaks
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
        /* This time we just create a {{{MutableMesh}}} and use that to specify the spatial locations of cells.*/
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();  //**Changed**//

        /* We create the same number of cells as the previous test.*/
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), TRANSIT);

        /* This time we create a {{{MeshBasedCellPopulation}}} as we are using a {{{MutableMesh}}}.*/
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells); //**Changed**//

        /* We create an {{{OffSimulation}}} object as before, all we change is the output directory.*/
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo3"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);

        /* We use a different {{{Force}}} which is suitable for mesh based simulations.
         * Note we could of used the same one as before as node based and mesh based simulations
         * share many of the same forces the only difference between the two models is in how cell-cell interactions
         * are specified.
         */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force); //**Changed**//
        simulator.AddForce(p_force);

        /* Again we call the {{{Solve}}} method on the simulation to run the simulation.*/
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
        /* This time we just create a {{{MutableMesh}}} and use that to specify the spatial locations of cells.
         * Here we pass an extra argument to the {{{HoneycombMeshGenerator}}} which adds another 2 rows of
         * nodes round the mesh, known as ghost nodes.*/
        HoneycombMeshGenerator generator(2, 2, 2); //**Changed**//
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* We only want to create cells for non ghost nodes. To find these we get them from the {{{HoneycombMeshGenerator}}}
         * using the method {{{GetCellLocationIndices}}}. We also use a different {{{CellCycleModel}}}. Here we use a
         * {{{TysonNovakCellCycleModel}}} which solves a coupled set of ODEs for each cell to calculate when each cell divides. */
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();//**Changed**//
        std::vector<CellPtr> cells;
        CellsGenerator<TysonNovakCellCycleModel, 2> cells_generator; //**Changed**//
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT); //**Changed**//

        /* This time we create a {{{MeshBasedCellPopulation}}} as we are using a {{{MutableMesh}}} and have ghost nodes.
         * We also need to pass the indices of non ghost nodes as an extra argument.*/
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); //**Changed**//

        /* We create an {{{OffSimulation}}} object as before, all we change is the output directory and the end time.
         * The Tyson Novak model is for yeast cells and therefore cells proliferate much more often and so we run the simulation for
         * less time to keep cell numbers relatively small for this demo.
         *
         */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo4"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(2.0); //**Changed**//

        /* We use the same {{{Force}}} as before and run the simulation in the same way.*/
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
        /* We now want to impose periodic boundaries on the domain. To do this we create a {{{Cylindrical2dMesh}}}
         * using a {{{CylindricalHoneycombMeshGenerator}}}.*/
        CylindricalHoneycombMeshGenerator generator(5, 2, 2); //**Changed**//
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh(); //**Changed**//

        /* Again we create one cell for each non ghost node. Note that we have changed back to using a {{{StochasticDurationCellCycleModel}}}.*/
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator; //**Changed**//
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT);

        /* We use the same {{{CellPopulation}}}, {{{CellBasedSimulation}}} (only changing the output directory and end time) and {{{Force}}} as before and run the simulation.*/
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo5"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0); //**Changed**//

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
     * We next show how to modify the previous test to include one
     * or more 'obstructions' within the domain.
     */
    void TestMeshBasedMonolayerPeriodicSolidBottomBoundary() throw (Exception)
    {
        /* We make the same {{{Mesh}}}, {{{Cells}}}, {{{CellPopulation}}},
         * {{{CellBasedSimulation}}} and forces as before, all we change is the output directory.*/
        CylindricalHoneycombMeshGenerator generator(5, 2, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), STEM);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo6"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(20.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);


        /* We now want to impose the condition y>0 on the cells. To do this we create a "shared pointer" to a {{{PlaneBoundaryCondition}}}.
         * Much like the {{{RandomCellKiller}}} earlier we pass arguments to the constructor (a point (0,0) on the plane (line in 2D) and an outward pointing normal to the plane (0,-1) ) using the {{{MAKE_PTR_ARGS}}} macro.*/
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        /* Finally we call the {{{Solve}}} method as in all other simulations.*/
        simulator.Solve();
    }
    /*
     * EMPTYLINE
     *
     * The results may be visualized using {{{Visualize2dCentreCells}}} as described in the
     * previous test, with the results directory changed from {{{CellBasedDemo5}}} to {{{CellBasedDemo6}}}.
     *
     * EMPTYLINE
     *
     * == Test 7 - basic Potts-based simulation ==
     *
     * EMPTYLINE
     *
     * In the final test we show how to modify the earlier tests (using off lattice models) to implement a 'Potts-based' simulation,
     * in which cells are represented by collections of sites on a fixed lattice.
     */
   void TestPottsBasedMonolayer() throw (Exception)
   {
       /* In common with the off lattice simulations we begin by creating a mesh. Here we use the {{{PottsMeshGenerator}}}
        * class to generate a {{{PottsMesh}}} each element in the mesh is a collection of lattice sites (represented by nodes at their centres).
        * All the connectivity between lattice sites is defined by the {{{PottsMeshGenerator}}},
        * and there are arguments to make the domains periodic.
        */
       PottsMeshGenerator<2> generator(20, 2, 4, 20, 2, 4); //**Changed**//
       PottsMesh<2>* p_mesh = generator.GetMesh(); //**Changed**//

       /* We generate one cell for each element as in vertex based simulations.*/
       std::vector<CellPtr> cells;
       CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
       cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), TRANSIT);

       /* As we have a {{{PottsMesh}}} we use a {{{PottsBasedCellPopulation}}}. Note here we also change the
        * "temperature" of the Potts simulation to make cells more motile.*/
       PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);//**Changed**//
       cell_population.SetTemperature(1.0);

       /* As a Potts simulation is restricted to a lattice we create a {{{OnSimulation}}} object and pass in the {{{CellPopulation}}} in much the same
        * way as an {{{OffLatticeSimulation}}} in the above examples. We also set some
        * options on the simulation like output directory and end time.
        */
       OnLatticeSimulation<2> simulator(cell_population);//**Changed**//
       simulator.SetOutputDirectory("CellBasedDemo7"); //**Changed**//
       simulator.SetEndTime(20.0);

       /* In order to specify how cells move around we create "shared pointer"s to
        * {{{UpdateRule}}} objects and pass them to the {{{OnLatticeSimulation}}}.
        * This is analogous to {{{Forces}}} in earlier examples.
        */
       MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule); //**Changed**//
       simulator.AddPottsUpdateRule(p_volume_constraint_update_rule); //**Changed**//
       MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule); //**Changed**//
       simulator.AddPottsUpdateRule(p_surface_area_update_rule); //**Changed**//
       MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule); //**Changed**//
       simulator.AddPottsUpdateRule(p_adhesion_update_rule); //**Changed**//

       /* We can add {{{CellKillers}}} as before.*/
       MAKE_PTR_ARGS(RandomCellKiller<2>, p_cell_killer, (&cell_population, 0.01));
       simulator.AddCellKiller(p_cell_killer);

       /* Again we run the simulation by calling the {{{Solve}}} method.*/
       simulator.Solve();
   }
   /*
    * EMPTYLINE
    *
    * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
    * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/CellBasedDemo7/results_from_time_0}}}.
    * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
    * java executable.
    *
    * EMPTYLINE
    */
};

#endif /*TESTCELLBASEDDEMOTUTORIAL_HPP_*/