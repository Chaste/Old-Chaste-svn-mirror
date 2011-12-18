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

#ifndef TESTRUNNINGPOTTSBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGPOTTSBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and visualize Potts-based simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize Potts-based simulations.
 * Full details of the mathematical model proposed by Graner, F. and Glazier, J. A. (1992). Simulation
 * of biological cell sorting using a two-dimensional extended potts model. Phys. Rev. Lett., 69(13):2013–2016.
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

/* The remaining header files define classes that will be used in the cell population
 * simulation test. We have encountered some of these header files in previous cell-based
 * Chaste tutorials. */
#include "CellsGenerator.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the cell-cycle model,one with stochastic cell-cycle times. */
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable mesh. */
#include "PottsMeshGenerator.hpp"
/* The next header file defines the class that simulates the evolution of an on lattice {{{CellPopulation}}}. */
#include "OnLatticeSimulation.hpp"
/* The next header file defines a potts-based {{{CellPopulation}}} class.*/
#include "PottsBasedCellPopulation.hpp"
/* The next header file defines some update rules for describing the Hamiltonain used to define the Potts simulations. */
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"

/*
 * Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}
 * and defines some test methods.
 */
class TestRunningPottsBasedSimulationsTutorial : public CxxTest::TestSuite
{
public:
    /* EMPTYLINE
    *
    * == Test 1 - a basic Potts-based simulation ==
    *
    * EMPTYLINE
    *
    * In the first test, we run a simple Potts-based simulation, in which we create a monolayer
    * of cells, using a Potts mesh. Each cell is assigned a fixed cell-cycle model.
    */

    void TestMonolayer() throw(Exception)
    {
        /* As in previous cell-based Chaste tutorials, we begin by setting up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Next, we generate a Potts mesh. To create a {{{PottsMesh}}}, we can use
        * the {{{PottsMeshGenerator}}}. This generates a regular square-shaped mesh,
        * in which all elements are the same size.
        *
        * Here the first three arguments specify the domain width; the number of elements across; and the width of
        * elements. The second set of three arguments specify the domain height; the number of elements up; and
        * the height of individual elements.
        * We have chosen a 2 by 2 block of elements each consisting of 4 by 4  ( = 16) lattice sites.
        */
        PottsMeshGenerator<2> generator(50, 2, 4, 50, 2, 4);  // Parameters are: lattice sites across; num elements across; element width; lattice sites up; num elements up; and element height
        PottsMesh<2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
        * To do this, we the `CellsGenerator` helper class, which is templated over the type
        * of cell model required (here {{{StochasticDurationGenerationBasedCellCycleModel}}})
        * and the dimension. We create an empty vector of cells and pass this into the
        * method along with the mesh. The second argument represents the size of that the vector
        * {{{cells}}} should become - one cell for each element. Third argument makes all cells
        * proliferate a specified numner of times (defaults to three).*/
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),TRANSIT);

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
        * In general, this class associates a collection of cells with a set of elements or a mesh.
        * For this test, because we have a {{{PottsMesh}}}, we use a particular type of
        * cell population called a {{{PottsBasedCellPopulation}}}.
        */
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /*
         * We can set the "Temperature" to be used in the Potts Simulation using the optional command.
         * The default value is 0.1.
         */
        cell_population.SetTemperature(0.1);
        /*
         * By default the Potts simulation will make 1 sweep over the whole domain per timestep.  To use a different
         * number of sweeps per timestep use the command.
         */
        cell_population.SetNumSweepsPerTimestep(1);


        /* We then pass in the cell population into a {{{OnLatticeSimulation}}},
         * and set the output directory and end time.*/
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("PottsBasedMonolayer");
        simulator.SetEndTime(50.0);
        /*
         * The default timestep is 0.1, but can be changed using the below command. The timestep is used in conjuction with the "Temperature" and
         * number of sweeps per timestep to specify the relationship between cell movement and prolideration.
         */
        simulator.SetDt(0.1);

        /* We must now create one or more update rules, which determine the Hamiltonian
        * in the Potts simulation. For this test, we use two update rules based upon
        * an area constraint and adhesion between cells and pass them to the {{{OnLatticeSimulation}}}.
        * For a list of possible update rules see subclasses of {{{AbstractPottsUpdateRule}}}.
        * These can be found in the inheritance diagram, here, [class:AbstractPottsUpdateRule AbstractPottsUpdateRule].
        *
        * Similarly to specifying forces for off lattice simulations we use the {{{MAKE_PTR}}} macro
        * to make a pointer to our required update rule before specifying parameters and passing to the simulation as follows
        */
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        /*
         * Set an appropriate target volume in number of lattice sites. Here we use the default vale of 16 lattice sites.
         */
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        /*
         * You can also vary the deformatuion energy parameter. The larger the parameter
         * the more cells will try to maintain target volume. Here we use the default vale of 0.2.
         */
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        /*
         * Finally we add the update rule to the simulator.
         */
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        /*
         * We repeat the process for any other update rules.
         */
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* {{{SimulationTime::Destroy()}}} '''must''' be called at the end of the test.
        * If not, when {{{SimulationTime::Instance()->SetStartTime(0.0);}}} is called
        * at the beginning of the next test in this file, an assertion will be triggered.
        */
        SimulationTime::Destroy();
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
    * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/PottsBasedMonolayer/results_from_time_0}}}.
    * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
    * java executable.
    *
    * EMPTYLINE
    *
    * == Test 2 - Cell Sorting ==
    *
    * EMPTYLINE
    *
    * The next test generates a collection of cells, there are two types of cells labelled ones
    * and non labelled ones, there is differential adhesion between the cell types. For the
    * specified parameters,the cells sort into separate types.
    *
    * Parameters are taken from Graner, F. and Glazier, J. A. (1992). Simulation of biological
    * cell sorting using a two-dimensional extended potts model. Phys. Rev. Lett., 69(13):2013–2016.
    *
    */

    void TestPottsMonolayerCellSorting() throw (Exception)
    {
        /* As in previous cell-based Chaste tutorials, we begin by setting up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Next, we generate a Potts mesh. To create a {{{PottsMesh}}}, we can use
        * the {{{PottsMeshGenerator}}}. This generates a regular square-shaped mesh,
        * in which all elements are the same size.
        *
        * Here the first three arguments specify the domain width; the number of elements across; and the width of
        * elements. The second set of three arguments specify the domain height; the number of elements up; and
        * the height of individual elements.
        * We have chosen an 8 by 8 block of elements each consisting of 4 by 4  ( = 16) lattice sites.
        */
        PottsMeshGenerator<2> generator(50, 8, 4, 50, 8, 4);  // Parameters are: lattice sites across; num elements across; element width; lattice sites up; num elements up; and element height
        PottsMesh<2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
        * To do this, we the `CellsGenerator` helper class, as before but this time the third argument is set to
        * DIFFERENTIATED to make all cells non-proliferative.*/
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        /* Before we make a {{{CellPopulation}}} we make a pointer to a cell label and then assign this
         * label to some randomly chosen cells. */
        MAKE_PTR(CellLabel, p_label);
        for (unsigned i = 0; i<cells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                cells[i]->AddCellProperty(p_label);
            }
        }

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
        * In general, this class associates a collection of cells with a set of elements or a mesh.
        * For this test, because we have a {{{PottsMesh}}}, we use a particular type of
        * cell population called a {{{PottsBasedCellPopulation}}}.
        */
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We then pass in the cell population into a {{{OnLatticeSimulation}}},
         * and set the output directory and end time. */
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("PottsMonolayerCellSorting");
        simulator.SetEndTime(20.0);

        /* We must now create one or more update rules, which determine the Hamiltonian
        * in the Potts simulation. For this test, we use two update rules based upon
        * an area constraint and differential adhesion between cells and pass them to the {{{OnLatticeSimulation}}}.
        */
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        /*
         * Set an appropriate target volume in number of lattice sites.
         */
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<2>, p_differential_adhesion_update_rule);
        /*
         * These parameters cause the cells to sort, for different values you can get different patterns.
         */
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.11);
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddPottsUpdateRule(p_differential_adhesion_update_rule);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* {{{SimulationTime::Destroy()}}} '''must''' be called at the end of the test.
        * If not, when {{{SimulationTime::Instance()->SetStartTime(0.0);}}} is called
        * at the beginning of the next test in this file, an assertion will be triggered.
        */
        SimulationTime::Destroy();
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
    * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/PottsMonolayerCellSorting/results_from_time_0}}}.
    * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
    * java executable.
    *
    * EMPTYLINE
    *
    * == Test 3 - 3D Cell Sorting ==
    *
    * EMPTYLINE
    *
    * The next test extends the previous example to three dimensions.
    *
    */

    void TestPottsSpheroidCellSorting() throw (Exception)
    {
        /* As in previous cell-based Chaste tutorials, we begin by setting up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Next, we generate a Potts mesh. To create a {{{PottsMesh}}}, we can use
        * the {{{PottsMeshGenerator}}}. This generates a regular square-shaped mesh,
        * in which all elements are the same size.
        *
        * Here the first three arguments specify the domain width; the number of elements across; and the width of
        * elements. The second set of three arguments specify the domain height; the number of elements up; and
        * the height of individual elements.  The third set of three arguments specify the domain depth; the number of elements deep; and
        * the depth of individual elements.
        * We have chosen an 4 by 4 by 4 ( = 32) block of elements each consisting of 2 by 2 by 2 ( = 8) lattice sites.
        */
        PottsMeshGenerator<3> generator(10, 4, 2, 10, 4, 2, 10, 4, 2);  // Parameters are: lattice sites across; num elements across; element width; lattice sites up; num elements up; and element height; lattice sites deep; num elements deep; and element depth
        PottsMesh<3>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
        * To do this, we the `CellsGenerator` helper class, as before but this time the third argument is set to
        * DIFFERENTIATED to make all cells non-proliferative.*/
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        /* As for the 2D case before we make a {{{CellPopulation}}} we make a pointer to a cell label and then assign this
         * label to some randomly chosen cells. */
        MAKE_PTR(CellLabel, p_label);
        for (unsigned i = 0; i<cells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                cells[i]->AddCellProperty(p_label);
            }
        }

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
        * In general, this class associates a collection of cells with a set of elements or a mesh.
        * For this test, because we have a {{{PottsMesh}}}, we use a particular type of
        * cell population called a {{{PottsBasedCellPopulation}}}.
        */
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);

        /* We then pass in the cell population into a {{{OnLatticeSimulation}}},
         * and set the output directory and end time. */
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("PottsCellSorting3D");
        simulator.SetEndTime(20.0);

        /* We must now create one or more update rules, which determine the Hamiltonian
        * in the Potts simulation. For this test, we use two update rules based upon
        * an area constraint and differential adhesion between cells and pass them to the {{{OnLatticeSimulation}}}.
        */
        MAKE_PTR(VolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        /*
         * Now set the target volume to be appropriate for this 3D simulation.
         */
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(8.0);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);

        /*
         * We use the same differential adhesion parameters as in the 2D case.
         */
        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<3>, p_differential_adhesion_update_rule);
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.11);
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddPottsUpdateRule(p_differential_adhesion_update_rule);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* {{{SimulationTime::Destroy()}}} '''must''' be called at the end of the test.
        * If not, when {{{SimulationTime::Instance()->SetStartTime(0.0);}}} is called
        * at the beginning of the next test in this file, an assertion will be triggered.
        */
        SimulationTime::Destroy();
   }
    /*
       * EMPTYLINE
       *
       * To visualize the results, we need to use paraview see the UserTutorials/VisualizingWithParaview tutorial for more information.
       *
       * In particular we need to create glyphs for each lattice site and color them based on cell Index or cell type.
       *
       *
       * EMPTYLINE
       */


};

#endif /* TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_ */
