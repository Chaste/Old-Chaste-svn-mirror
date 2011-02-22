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
 * a vertex model. To be able to view these simulations, you must first have
 * downloaded and installed VTK and Paraview, and updated your hostconfig file
 * to ensure that it knows to use VTK.
 *
 * For the tests we require the following headers. Firstly, we need the test suite below,
 * which allows us to use certain methods in our test (this header file should be included
 * in any Chaste test).
 */
#include <cxxtest/TestSuite.h>

/* Any test in which the {{{GetIdentifier()}}} method is used, 
 * even via the main cell_based code ({{{AbstractCellPopulation}}} output methods), must 
 * include {{{CheckpointArchiveTypes.hpp}}} 
 * or {{{CellBasedSimulationArchiver.hpp}}} as the first Chaste header included. 
 */
#include "CheckpointArchiveTypes.hpp" 


/*
 * The next two header files define a stochastic and fixed duration cell-cycle model respectively.
 */
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable mesh for a cell-centre model. */
#include "HoneycombMeshGenerator.hpp"
/* The next header file defines a helper class for generating a suitable vertex mesh. */
#include "HoneycombVertexMeshGenerator.hpp"
/* The next header file defines a helper class for generating
 * a vector of cells for a given mesh. */
#include "CellsGenerator.hpp"
/* The next header file defines a {{{CellPopulation}}} class that uses a mesh, which contains ghost nodes. */
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
/* The next header file defines a vertex-based {{{CellPopulation}}} class.*/
#include "VertexBasedCellPopulation.hpp"
/* The next header file defines a force law, based on a linear spring, for describing
 * the mechanical interactions between neighbouring cells in the cell population.
 */
#include "GeneralisedLinearSpringForce.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population, subject to each vertex.
 */
#include "NagaiHondaForce.hpp"
/* The next header file defines the class that simulates the evolution of a {{{CellPopulation}}}.
 */
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
        /* As in '''all''' cell-based simulations, we must first set the start time.
         */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Next, we generate a mesh we use the {{{HoneycombMeshGenerator}}}. This
         * generates a honeycomb-shaped mesh, in which all nodes are equidistant.
         * Here the first and second arguments define the size of the mesh - we have
         * chosen a mesh that is 10 nodes (i.e. cells) wide, and 10 nodes high.
         * The third argument indicates that we want 2 layers of ghost nodes around the mesh.
         * The last boolean parameter indicates that we do not want cylindrical boundary
         * conditions. We generate a pointer to the mesh, and then get the location indices of
         * the real cells.
         */
        HoneycombMeshGenerator generator(10, 10, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * Then we loop over the number of real nodes in the mesh and assign a cell
         * to each node. Each cell will have a randomly chosen birth time. */
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(TRANSIT);
            p_model->SetMaxTransitGenerations(UINT_MAX);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (p_model->GetStemCellG1Duration()
                                    + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Now we have a mesh and a set of cells to go with it we can create a ''CellPopulation''.
         * In general, this class associates a collection of cells with a set of nodes or a mesh.
         * For this test we use a particular type of cell population called a
         * {{{MeshBasedCellPopulationWithGhostNodes}}}.
         */
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /* In order to output the .vtu files required for Paraview, we explicitly
         * instruct the simulation to output the data we need.
         */
        cell_population.SetOutputVoronoiData(true);

        /* Now we define the cell-based simulation object, passing in the cell population. */
        CellBasedSimulation<2> simulator(cell_population);

        /* Set the output directory on the simulator (relative to
         * "/tmp/<USER_NAME>/testoutput") and the end time (in hours).
         */
        simulator.SetOutputDirectory("Test2DMonolayerSimulationForVisualizing");
        simulator.SetEndTime(1.0);

        /* We must now create one or more force laws, which determine the mechanics of
         * the cell population. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring. We put a pointer
         * to this force into a vector. We use a cut-off point which represents that cells farther
         * than 1.5 cell lengths apart, do not exert forces on one another.
         */

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
	* To visualize the results, you must first open Paraview. Open the folder containing your test output using the 'file' menu at
	* the top. The output will be located in {{{/tmp/$USER/testoutput/Test2DMonolayerSimulationForVisualizing/results_from_time_0}}}.
	* There will be a .vtu file generated for every timestep, which must all be opened at once to view the simulation. To do this,
	* simply select {{{results_..vtu}}}. You should now see {{{results_*}}} in the pipeline browser. Click {{{Apply}}} in the properties tab
	* of the object inspector, and you should now see a visualisation in the right hand window.
	*
	* At this stage, it will be necessary to refine how you wish to view this particular visualisation. The viewing styles can be edited using
	* the display tab of the object inspector. In particular, under {{{Style}}}, the representation drop down menu allows you to view
	* the cells as a surface with edges, or as simply a wireframe. It is advisable at this point to make yourself familiar with the different
	* viewing options, colour and size settings.
	*
	* At this stage, the viewer is showing all cells in the simulation, including the ghost nodes. In order to view only real cells, you must
	* apply a threshold. This is achieved using the threshold button on the third toolbar (the icon is a cube with a green 'T' inside). Once you
	* click the threshold button, you will see a new threshold appear below your results in the pipeline browser. Go to the properties tab and
	* reset the lower threshold to be less than 0, and the upper threshold to be between 0 and 1, ensuring that the 'Non-ghosts' option is
	* selected in the 'Scalars' drop down menu. Once you have edited this, click apply (you may need to click it twice), and the visualisation on the
	* right window will have changed to eliminate ghost nodes.
	*
	* To view the simulation, simply use the animation buttons located on the top toolbar. You can also save a screenshot, or an animation, using
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

		/* Next, we generate a vertex mesh. To create a {{{MutableVertexMesh}}}, we can use
		* the {{{HoneycombVertexMeshGenerator}}}. This generates a honeycomb-shaped mesh,
		* in which all nodes are equidistant. Here the first and second arguments
		* define the size of the mesh - we have chosen a mesh that is 6 elements (i.e.
		* cells) wide, and 9 elements high.
		*/
		HoneycombVertexMeshGenerator generator(6, 9);	// Parameters are: cells across, cells up
		MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

		/* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
		* To do this, we the `CellsGenerator` helper class, which is templated over the type
		* of cell model required (here {{{FixedDurationGenerationBasedCellCycleModel}}})
		* and the dimension. We create an empty vector of cells and pass this into the
		* method along with the mesh. The second argument represents the size of that the vector
		* {{{cells}}} should become - one cell for each element. */
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

		/* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
		* In general, this class associates a collection of cells with a set of elements or a mesh.
		* For this test, because we have a {{{MutableVertexMesh}}}, we use a particular type of
		* cell population called a {{{VertexBasedCellPopulation}}}.
		*/
		VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

		/* Now we define the cell-based simulation object, passing in the cell population and collection
		* of force laws:
		*/
		CellBasedSimulation<2> simulator(cell_population);

		/* Set the output directory on the simulator and the end time (in hours).
		*/
		simulator.SetOutputDirectory("Test2DVertexMonolayerSimulationForVisualizing");
		simulator.SetEndTime(1.0);

        /* We must now create one or more force laws, which determine the mechanics of the vertices
         * of each cell in a cell population. For this test, we use one force law, based on the
         * Nagai-Honda mechanics. We put a pointer to this force into a vector.
         */

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
	* To visualize the results, follow the instructions above for the first simulation, ensuring that you open the
	* test output from the new folder, {{{Test2DVertexMonolayerSimulationForVisualizing}}}.
	*
	*/

};


#endif /* TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_ */
