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
#ifndef TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_
#define TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_
/*
 * = Examples showing how to visualize simulations in Paraview =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste is used to generate simulations
 * that can be viewed in Paraview, and how to use Paraview itself. Two examples
 * are provided - one using a cell-centre based model, and the second using
 * a vertex model. To be able to view these simulations, we must first have
 * downloaded and installed VTK and Paraview, and updated our hostconfig file
 * to ensure that it knows to use VTK.
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
 * simulation test. We have encountered each of these header files in previous cell-based
 * Chaste tutorials. */
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "CellBasedSimulation.hpp"

/* Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}
 * and defines some test methods.
 */
class TestVisualizingWithParaviewTutorial : public CxxTest::TestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1 - a cell-based monolayer simulation ==
     *
     * EMPTYLINE
     *
     * In the first test, we run a simple cell-based simulation, in which we use
     * a honeycomb mesh with ghost nodes, and give each cell a stochastic cell-cycle model.
     */
    void Test2DMonolayerSimulationForVisualizing() throw (Exception)
    {
        /* As in previous cell-based Chaste tutorials, we begin by setting up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* In a similar way to previous cell-based Chaste tutorials,
         * we create a mesh-based cell population in which cells are defined by their centres,
         * and cell proliferation is governed by a stochastic generation-based cell-cycle model
         * with no differentiation.
         */
        HoneycombMeshGenerator generator(10, 10, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(TRANSIT);
            p_model->SetMaxTransitGenerations(UINT_MAX);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()* (p_model->GetStemCellG1Duration() + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /*
         * The following line tells the cell population to write data to .vtu files with cells
         * not as points, but as polytopes. This is the default setting: we include the call
         * here to highlight this option. If writing point data, we may choose the shape used
         * to visualize each cell in Paraview using glyphs.
         */
        cell_population.SetWriteVtkAsPoints(false);

        /* In order to output the .vtu files required for Paraview, we explicitly
         * instruct the simulation to output the data we need.
         */
        cell_population.SetOutputVoronoiData(true);

        /* We then pass in the cell population into a {{{CellBasedSimulation}}},
         * and set the output directory and end time. */
        CellBasedSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test2DMonolayerSimulationForVisualizing");
        simulator.SetEndTime(1.0);

        /* We create a force law and pass it to the {{{CellBasedSimulation}}}. */
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* {{{SimulationTime::Destroy()}}} and {{{RandomNumberGenerator::Destroy()}}} '''must''' be called at the end of the test.
         * If not, when {{{SimulationTime::Instance()->SetStartTime(0.0);}}} is called
         * at the beginning of the next test in this file, an assertion will be triggered.
         */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, we must first open Paraview. We open the folder containing our test output using the 'file' menu at
    * the top. The output will be located in {{{/tmp/$USER/testoutput/Test2DMonolayerSimulationForVisualizing/results_from_time_0}}}.
    * There will be a .vtu file generated for every timestep, which must all be opened at once to view the simulation. To do this,
    * simply select {{{results_..vtu}}}. We should now see {{{results_*}}} in the pipeline browser. We click {{{Apply}}} in the properties tab
    * of the object inspector, and we should now see a visualization in the right hand window.
    *
    * At this stage, it will be necessary to refine how we wish to view this particular visualisation. The viewing styles can be edited using
    * the display tab of the object inspector. In particular, under {{{Style}}}, the representation drop down menu allows us to view
    * the cells as a surface with edges, or as simply a wireframe. It is advisable at this point to familiarize ourselves with the different
    * viewing options, colour and size settings.
    *
    * At this stage, the viewer is showing all cells in the simulation, including the ghost nodes. In order to view only real cells, we must
    * apply a threshold. This is achieved using the threshold button on the third toolbar (the icon is a cube with a green 'T' inside). Once you
    * click the threshold button, you will see a new threshold appear below your results in the pipeline browser. Go to the properties tab and
    * reset the lower threshold to be less than 0, and the upper threshold to be between 0 and 1, ensuring that the 'Non-ghosts' option is
    * selected in the 'Scalars' drop down menu. Once we have edited this, we click apply (we may need to click it twice), and the visualisation on the
    * right window will have changed to eliminate ghost nodes.
    *
    * To view the simulation, simply use the animation buttons located on the top toolbar. We can also save a screenshot, or an animation, using
    * the appropriate options from the file menu. Next to the threshold button are two other useful options, 'slice' and 'clip', but these will
    * only be applicable for 3D visualisations.
    *
    * EMPTYLINE
    *
    * == Test 2 - a basic vertex-based simulation ==
    *
    * EMPTYLINE
    *
    * Here, we run a simple vertex-based simulation, in which we create a monolayer
    * of cells using a mutable vertex mesh. Each cell is assigned a fixed cell-cycle model.
    */
    void TestMonolayerFixedCellCycle() throw(Exception)
    {
        /* First re-initialize time to zero. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* In this test, we create a vertex-based cell population in which cells are defined
         * by their vertices, and cell proliferation is governed by a fixed generation-based
         * cell-cycle model (with differentiation after a default number of generations).
         */
        HoneycombVertexMeshGenerator generator(6, 9);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We then pass in the cell population into a {{{CellBasedSimulation}}},
         * and set the output directory and end time. */
        CellBasedSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test2DVertexMonolayerSimulationForVisualizing");
        simulator.SetEndTime(1.0);

        /* We create a force law and pass it to the {{{CellBasedSimulation}}}. */
        NagaiHondaForce<2> nagai_honda_force;
        simulator.AddForce(&nagai_honda_force);

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
    * To visualize the results, we follow the instructions above for the first simulation, ensuring that we open the
    * test output from the new folder, {{{Test2DVertexMonolayerSimulationForVisualizing}}}.
    *
    */
};

#endif /* TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_ */
