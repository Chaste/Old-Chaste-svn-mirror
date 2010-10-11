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
/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and visualize vertex-based simulations on periodic meshes with different cell cycle models =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste is used to create, run and visualize vertex-based simulations.
 * Full details of the mechanical model proposed by T. Nagai and H. Honda ("A dynamic cell model for
 * the formation of epithelial tissues", Philosophical Magazine Part B 81:699-719).
 *
 * The first thing to do is include the following header, which allows us
 * to use certain methods in our test (this header file should be included
 * in any Chaste test).
 */
#include <cxxtest/TestSuite.h>
/* Any test in which the {{{GetIdentifier()}}} method is used, 
 * even via the main cell_based code ({{{AbstraceCellPopulation}}} output methods), must 
 * include {{{CheckpointArchiveTypes.hpp}}} 
 * or {{{CellBasedSimulationArchiver.hpp}}} as the first Chaste header included. 
 */
#include "CheckpointArchiveTypes.hpp" 


/* The next header file defines a helper class for generating
 * a vector of cells for a given mesh. */
#include "CellsGenerator.hpp"
/* The next header file defines a helper class for generating
 * cells for crypt simulations. */
#include "CryptCellsGenerator.hpp"
/*
 * The next three header files define three different types of cell-cycle model,
 * one with fixed cell-cycle times, one with stochastic cell-cycle times and one
 * where the cell-cycle time depends on the Wnt concentration.
 */
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable mesh. */
#include "HoneycombMutableVertexMeshGenerator.hpp"
/* The next header file defines a helper class for generating a periodic vertex mesh. */
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
/* The next header file defines a Wnt singleton class, which (if used) deals with the
 * imposed Wnt gradient in our crypt model.
 */
#include "WntConcentration.hpp"
/* The next header file defines the class that simulates the evolution of a {{{CellPopulation}}} */
#include "CellBasedSimulation.hpp"
/* The next header file defines the class that simulates the evolution of a crypt {{{CellPopulation}}}
 * for a vertex mesh. */
#include "VertexCryptSimulation2d.hpp"
/* The next header file defines a vertex-based {{{CellPopulation}}} class.*/
#include "VertexBasedCellPopulation.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population, subject to each vertex.
 */
#include "NagaiHondaForce.hpp"
/* The final header file defines a cell killer class, which implements sloughing of cells
 * into the lumen once they reach the top of the crypt.
 */
#include "SloughingCellKiller.hpp"
/* Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}
 * and defines some test methods.
 */
class TestRunningVertexBasedSimulationsTutorial : public CxxTest::TestSuite
{
public:
	/* EMPTYLINE
	*
	* == Test 1 - a basic vertex-based simulation ==
	*
	* EMPTYLINE
	*
	* In the first test, we run a simple vertex-based simulation, in which we create a monolayer
	* of cells, using a mutable vertex mesh. Each cell is assigned a fixed cell-cycle model.
	*/
	void TestMonolayerFixedCellCycle() throw(Exception)
	{
    	/* As in '''all''' cell-based simulations, we must first set the start time.
    	*/
    	SimulationTime::Instance()->SetStartTime(0.0);
    
    	/* Next, we generate a vertex mesh. To create a {{{MutableVertexMesh}}}, we can use
    	* the {{{HoneycombMutableVertexMeshGenerator}}}. This generates a honeycomb-shaped mesh,
    	* in which all nodes are equidistant. Here the first and second arguments
    	* define the size of the mesh - we have chosen a mesh that is 6 elements (i.e.
    	* cells) wide, and 9 elements high.
    	*/
    	HoneycombMutableVertexMeshGenerator generator(6, 9);	// Parameters are: cells across, cells up
    	MutableVertexMesh<2,2>* p_mesh = generator.GetMutableMesh();
    
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
    
    	/* We must now create one or more force laws, which determine the mechanics of the vertices
    	* of each cell in a cell population. For this test, we use one force law, based on the
    	* Nagai-Honda mechanics. We put a pointer to this force into a vector.
    	*/
    	NagaiHondaForce<2> force;
    	std::vector<AbstractForce<2>* > force_collection;
    	force_collection.push_back(&force);
    
    	/* Now we define the cell-based simulation object, passing in the cell population and collection
    	* of force laws:
    	*/
    	CellBasedSimulation<2> simulator(cell_population, force_collection);
    
    	/* Set the output directory on the simulator (relative to
    	* "/tmp/<USER_NAME>/testoutput") and the end time (in hours).
    	*/
    	simulator.SetOutputDirectory("MonolayerFixedCellCycle");
    	simulator.SetEndTime(1.0);
    
    	/* For longer simulations, you may not want to output the results
    	* every time step. In this case you can use the following method,
    	* to print results every 10 time steps instead. As the time step
    	* used by the simulator, is 30 s, this method will cause the
    	* simulator to print results every 5 min.
    	*/
    	//simulator.SetSamplingTimestepMultiple(10);
    
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
	* then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/MonolayerFixedCellCycle/results_from_time_0}}}.
	* You may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
	* java executable.
	*
	* EMPTYLINE
	*
	* When you visualize the results, you should see the cells whose centres lie at and above 4.0 dividing first. This is due
	* to the implementation of the {{{CellsGenerator}}}, which assigned a birthtime of (0 - i), where i is the element index of the cell.
	*
	* EMPTYLINE
	*
	* == Test 2 - create a vertex-based crypt simulation ==
	*
	* EMPTYLINE
	*
	* The next test generates a crypt, in which we use a cylindrical vertex mesh,
	* give each cell a fixed cell-cycle model, and enforce sloughing at the top of
	* the crypt.
	*/
	void TestVertexBasedCrypt() throw(Exception)
	{
	    /* First re-initialize time to zero. */
	    SimulationTime::Instance()->SetStartTime(0.0);

	    /* Create a cylindrical mesh, and get the cell location indices. To enforce
	     * periodicity at the left and right hand sides of the mesh, we use a subclass
	     * called {{{Cylindrical2dMesh}}}, which has extra methods for maintaining
	     * periodicity.
	     */
  	    CylindricalHoneycombVertexMeshGenerator generator(6, 9);
 	    Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
    	 * To do this, we the `CryptCellsGenerator` helper class, which is templated over the type
    	 * of cell model required (here {{{FixedDurationGenerationBasedCellCycleModel}}})
    	 * and the dimension. We create an empty vector of cells and pass this into the
    	 * method along with the mesh. The third argument 'true' indicates that the cells
    	 * should be assigned random birth times, to avoid synchronous division. The
    	 * {{{cells}}} vector is populated once the method {{{Generate}}} is
    	 * called.
	     * The last four arguments represent the height below which cells belong to generations 0,
	     * 1, 2, 3 and 4, respectively.
	     */
	    std::vector<CellPtr> cells;
	    CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
	    cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 1.0, 2.0, 3.0, 4.0);

    	/* Create cell population, as before. */
    	VertexBasedCellPopulation<2> crypt(*p_mesh, cells);
    
    	/* Create force law and force collection, as above. */
    	NagaiHondaForce<2> force_law;
	    std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

	    /* Create a simulator as before (except setting a different output directory). */
	    VertexCryptSimulation2d simulator(crypt, force_collection);
	    simulator.SetOutputDirectory("VertexCrypt");
        simulator.SetEndTime(1);

        /* Before running the simulation, we add a cell killer. This object
	     * dictates conditions under which cells die. For this test, we use
	     * a {{{SloughingCellKiller}}}, which kills cells above a certain height.
	     */
    	double crypt_length = 6.0;
    	SloughingCellKiller<2> sloughing_cell_killer(&crypt, crypt_length);
    	simulator.AddCellKiller(&sloughing_cell_killer);

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
	* then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/VertexCrypt/results_from_time_0}}}.
	* You may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
	* java executable.
	*
	* EMPTYLINE
	*
	* When you visualize the results, you should see three colours of cells: a row of blue stem cells, 3 rows of yellow transit
	* cells, and 5 rows of pink differentiated cells. Cells above 6.0 will be sloughed off immediately.
	*
	* EMPTYLINE
	*
	* == Test 3 - create a vertex-based crypt simulation with a simple wnt dependent cell cycle model ==
	*
	* EMPTYLINE
	*
	* The next test generates a crypt, in which we use a cylindrical vertex mesh, and
	* impose a linearly decreasing concentration gradient of Wnt. Cells detect the level of Wnt
	* at their centre and those that are in a region of sufficient Wnt are defined to be transit cells,
	* whilst those above this Wnt threshold are defined to be differentiated. The cell cycle length of
	* transit cells is then assigned randomly from a uniform distribution.
	*/
	void TestVertexBasedCryptWithSimpleWntCellCycleModel() throw(Exception)
	{
	/* First re-initialize time to zero. */
	SimulationTime::Instance()->SetStartTime(0.0);

	/* Create a cylindrical mesh, and get the cell location indices, as before.
	*/
	CylindricalHoneycombVertexMeshGenerator generator(6, 9);
	Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

	/* Create a {{{std::vector}}} of {{{CellPtr}}}s.
	* Generate cells, which are assigned a {{{SimpleWntCellCycleModel}}} using
	* the {{{CryptCellsGenerator}}}. The final boolean argument 'true' indicates
	* to assign randomly chosen birth times.
	*/
	std::vector<CellPtr> cells;
	CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
	cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

	/* Create cell population, as before. */
	VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

    /* Set the crypt length this will be used for sloughing and calculating the Wnt gradient */
    double crypt_length = 6.0;

	/* The other change needed: Cells with a Wnt-based cell cycle need to know
	* the concentration of Wnt wherever they are. To do this, we set up a {{{WntConcentration}}}
	* class. This is another singleton class (ie accessible from anywhere), so all
	* cells and cell cycle models can access it. We need to say what the profile of the
	* Wnt concentation should be - here, we say it is {{{LINEAR}}} (linear decreasing from 1 to 0
	* from the bottom of the crypt to the top). We also need to inform the {{{WntConcentration}}}
	* of the cell population.*/
	WntConcentration<2>::Instance()->SetType(LINEAR);
	WntConcentration<2>::Instance()->SetCellPopulation(crypt);
	WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

	/* Create force law and force collection, as above. */
	NagaiHondaForce<2> force_law;
	std::vector<AbstractForce<2>*> force_collection;
	force_collection.push_back(&force_law);

	/* Create a simulator as before (except setting a different output directory). */
	VertexCryptSimulation2d simulator(crypt, force_collection);
	simulator.SetOutputDirectory("VertexCryptWithSimpleWntCellCycleModel");
	simulator.SetEndTime(1);

	/* Before running the simulation, we add a cell killer, as before.*/
	SloughingCellKiller<2> sloughing_cell_killer(&crypt, crypt_length);
	simulator.AddCellKiller(&sloughing_cell_killer);

	/* Here we impose a boundary condition at the base: that cells
	* at the bottom of the crypt are repelled if they move past 0.*/
	simulator.UseJiggledBottomCells();

	/* Run the simulation, by calling {{{Solve()}}}. */
	simulator.Solve();

	/* {{{SimulationTime::Destroy()}}} '''must''' be called at the end of the test.*/
	 SimulationTime::Destroy();
	}
	/*
	* EMPTYLINE
	*
	* To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
	* then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/VertexCryptWithSimpleWntCellCycleModel/results_from_time_0}}}.
	* You may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
	* java executable.
	*
	* EMPTYLINE
	*
	* When you visualize the results, you should see two colours of cells: yellow transit
	* cells and pink differentiated cells. Cells above 6.0 will be sloughed off immediately.
	*/

};

#endif /* TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_ */
