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

#ifndef TESTRUNNIGCONTACTINHIBITIONSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNIGCONTACTINHIBITIONSIMULATIONSTUTORIAL_HPP_

/*
 * = An example showing how to use the contact inhibition cell cycle model (with the contact inhibition simulator) =
 *
 * == Introduction ==
 *
 * In this tutorial, we will show how to use a simple implementation of the contact inhibition cell-cycle mode,
 * that stops cell division when the volume of the cell is smaller than a critical value.
 *
 * Firstly, we consider two mesh-based populations in 2-D with cells trapped in a square box. In the first population,
 * all the cells are contact inhibited and we study the effect of the critical volume upon the global cell density. In the
 * second population, we consider a mix of normal cells (contact inhibited) and tumour cells (not inhibited) and study the growth of the
 * tumour cells within the box.
 *
 * Secondly, we look at the behaviour of a vertex-based population in a box and the effect of contact inhibition.
 *
 * == Including header files ==
 *
 * We begin by including the necessary header files. */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* These two headers need to be includes here to ensure archiving of {{{CelwiseData}}} works on all Boost versions*/
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

/* The next header includes the Boost shared_ptr smart pointer, and defines some useful
 * macros to save typing when using it. */
#include "SmartPointers.hpp"
/* The next header include the NEVER_REACHED macro, used in one of the methods below. */
#include "Exception.hpp"

/*
 * The next header file defines the contact inhibition cell-cycle model that inherits from {{{AbstractSimpleCellCycleModel}}}.
 * The duration of the G1 phase depends on the deviation from a target volume (or area/length in 2D/1D): if the volume is
 * lower than a given fraction of the target volume, the G1 phase continues. The target volume and the critical fraction
 * are indicated in the user's Test file, and compared to the real volumes stored in {{{CellwiseData}}}, a singleton class.
 * This model allows for quiescence imposed by transient periods of high stress, followed by relaxation. Note that
 * in this cell cycle model, quiescence is implemented only by extending the G1 phase. Therefore, if a cell
 * is compressed during G2 or S phases then it will still divide, and thus cells whose volumes are smaller
 * than the given threshold may still divide.
 */
#include "ContactInhibitionCellCycleModel.hpp"

/*
 * The next header is the simulation class corresponding to the contact inhibition cell-cycle model. 
 * The essential difference with other simulators is that {{{CellWiseData}}} is updated with the 
 * volumes of the Voronoi elements representing each cell.
 */
#include "VolumeTrackedOffLatticeSimulation.hpp"
/* The remaining header files define classes that will be also be used and are presented in other tutorials. */
#include "MeshBasedCellPopulation.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "SimulationTime.hpp"
#include "CellLabel.hpp"
#include "MutableMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "PlaneBoundaryCondition.hpp"

/* We first define the global test class that inherits from {{{AbstractCellBasedTestSuite}}}. */
class TestRunningContactInhibitionSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    /*
     * == Testing healthy cell contact inhibition with mesh-based population ==
     *
     * In this first test we show how to simulate the behaviour of cells healthy cells trapped in a box.
     * Each cell will only divide if there is sufficient room.
     */
    void TestContactInhibitionInBox()
    {
        /* We use the honeycomb mesh generator to create a honeycomb mesh and
         * the associated mutable mesh. */
        HoneycombMeshGenerator generator(3, 3);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* We now create a vector of cell pointers. */
        std::vector<CellPtr> cells;

        /* We then define the mutation state of the cells we are working with. We will just consider
         * wild type mutations here. */
        MAKE_PTR(WildTypeCellMutationState, p_state);

        /* We now create a cell-cycle (only contact inhibited) model for these cells and loop over the
         * nodes of the mesh to create as many elements in the vector of cell pointers as there are
         * in the initial mesh. */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetCellProliferativeType(TRANSIT);
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-2.0*(double)i);
            p_cycle_model->SetQuiescentVolumeFraction(0.5);
            p_cycle_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        /* We now create a cell population, that takes several inputs: the mesh (for the position); and
         * the vector of cell pointers (for cycles and states)*/
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* In order to visualise labelled cells (i.e those that are inhibited from division) you need to use the following command.*/
        cell_population.SetOutputCellMutationStates(true);

        /* To keep track of the volumes of the cells that are used in the contact inhibition cell-cycle,
         * we use the singleton class {{{CellWiseData}}}. Here, we just initialise it with one variable
         * and associate it with the cell population. */

        /* This creates the instance of {{{CellwiseData}}}: */
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        /* the first thing we do is set the number of variables we with to use {{{CellWiseData}}} to track, we do this pay passing the number
         * of cells in the simulation and the number of variables to  the {{{SetNumCellsAndVars}}} method.*/
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        /* We then pass the cell population to the {{{CellWiseData}}} object so the size of the data can vary over time.
         * The simulation won't run without the cell population set.*/
        p_data->SetCellPopulation(&cell_population);
        /* We now loop over the cells and assign an initial value for the volume. This can be anything as it is overwritten
         * by the {{{PostSolve}}} method in {{{VolumeTrackedOffLatticeSimulation}}}.*/
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(1.0, cell_population.GetLocationIndexUsingCell(*cell_iter));
        }

        /* Then, we define the contact {{{VolumeTrackedOffLatticeSimulation}}} class, that automatically updates the volumes of the cells
         * in {{{CellWiseData}}}. We also set up the output directory, the end time and the output multiple.
         */
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestContactInhibitionInBox");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);

        /* Next, we create a force law (springs) to be applied between cell centres and set up a
         * cut-off length beyond which cells stop interacting. We then pass this to the {{{VolumeTrackedOffLatticeSimulation}}} */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        /*
         * To study the behaviour of the cells with varying volume, we trap them in a box, i.e., between
         *  4 plane boundary conditions. These planes are indicated by a point and a normal and then passed
         *  to the {{{VolumeTrackedOffLatticeSimulation}}}. The boundaries chosen are to make the test run
         *  in a short amount of time, if you can make the box larger then the test will take longer to run.
         */

        /* First x>0 */
        c_vector<double,2> point = zero_vector<double>(2);
		c_vector<double,2> normal = zero_vector<double>(2);
		point(0) = 0.0;
		point(1) = 0.0;
		normal(0) = -1.0;
		normal(1) = 0.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc1);
		/* Second x<2.5 */
		point(0) = 2.5;
		normal(0) = 1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc2);
		/* Third y>0 */
		point(0) = 0.0;
		point(1) = 0.0;
		normal(0) = 0.0;
		normal(1) = -1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc3);
		/* Finally y<2.5 */
		point(1) = 2.5;
		normal(1) = 1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc4);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* Finally, as in previous cell-based Chaste tutorials, we call {{{Destroy()}}} on the singleton classes. */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellwiseData<2>::Destroy();
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/TestContactInhibitionInBox/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * You will notice that once the cells are below a certain size they no longer proliferate and turn dark blue in the visualisation.
     *
     * EMPTYLINE
     *
     * == Testing normal and tumour cells with mesh-based population ==
     *
     * We now test the behaviour of a mixture of healthy and tumour cells in a Box. In this test healthy cells will only
     * divide if there is sufficient room whereas tumour cells will divide regardless.
     */
    void TestContactInhibitionInBoxWithMutants()
    {
        /* Create a simple mesh. */
        HoneycombMeshGenerator generator(3, 3);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* Create cell state. */
        MAKE_PTR(WildTypeCellMutationState, p_state);
        std::vector<CellPtr> cells;

        /* Create cells. The difference here is that one of the cells is not contact-inhibited, but rather
         * is defined by a {{{StochasticDurationCellCycleModel}}}. */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            if (i==1)
            {
                StochasticDurationCellCycleModel* p_cycle_model = new StochasticDurationCellCycleModel();
                p_cycle_model->SetCellProliferativeType(STEM);
                p_cycle_model->SetBirthTime(-14.0);

                CellPtr p_cell(new Cell(p_state, p_cycle_model));
                p_cell->SetBirthTime(0.0);
                cells.push_back(p_cell);
            }
            else
            {
                ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
                p_cycle_model->SetCellProliferativeType(TRANSIT);
                p_cycle_model->SetDimension(2);
                p_cycle_model->SetBirthTime(-2.0*(double)i);
                p_cycle_model->SetQuiescentVolumeFraction(0.8);
                p_cycle_model->SetEquilibriumVolume(1.0);

                CellPtr p_cell(new Cell(p_state, p_cycle_model));
                p_cell->InitialiseCellCycleModel();
                cells.push_back(p_cell);
            }
        }

        /* We now create a cell population, that takes several inputs: the mesh; and
         * the vector of cell pointers*/
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* In order to visualise labelled cells (i.e those that are inhibited from division) you need to use the following command.*/
        cell_population.SetOutputCellMutationStates(true);

        /* To keep track of the volumes of the cells that are used in the contact inhibition cell-cycle,
         * we use the singleton class {{{CellWiseData}}}. Here, we just initialise it with one variable
         * and associate it with the cell population: */
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(1.0, cell_population.GetLocationIndexUsingCell(*cell_iter));
        }

        /*  Then, we define the contact {{{VolumeTrackedOffLatticeSimulation}}} class, that automatically updates the volumes of the cells
         * in {{{CellWiseData}}}. We also set up the output directory, the end time and the output multiple.
         */
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestContactInhibitionTumourInBox");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);


        /* Next, we create a force law (springs) to be applied between cell centres and set up a
         * cut-off length beyond which cells stop interacting. We then pass this to the {{{VolumeTrackedOffLatticeSimulation}}} */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        /*
         *  Again we trap cells in a box, i.e., between 4 plane boundary conditions.
         *  These planes are indicated by a point and a normal and then passed
         *  to the {{{VolumeTrackedOffLatticeSimulation}}}. The boundaries chosen are to make the test run
         *  in a short amount of time, if you can make the box larger then the test will take longer to run.
         */

        /* First x>0 */
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = -1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        /* Second x<2.5 */
        point(0) = 2.5;
        normal(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);
        /* Third y>0 */
        point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        /* Finally y<2.5 */
        point(1) = 2.5;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* Finally, as in previous cell-based Chaste tutorials, we call {{{Destroy()}}} on the singleton classes. */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellwiseData<2>::Destroy();
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/TestContactInhibitionTumourInBox/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * You will notice that once the healthy cells (yellow) are below a certain size they no longer proliferate and turn dark blue in the visualisation.
     * Whereas Tumour cells (light blue) on the other hand will continue to proliferate.
     *
     * EMPTY LINE
     *
     * == Testing contact inhibition in vertex-based monolayer ==
	 *
	 * We now test the behaviour of normal contact inhibited cells for a vertex-based population.
	 * The example we use is a growing monolayer.
	 *
	 */
	void TestContactInhibitionWithVertex()
	{
		// Create a simple 2D MutableVertexMesh.
		HoneycombVertexMeshGenerator generator(2, 2);
		MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

		// Create cell state.
		MAKE_PTR(WildTypeCellMutationState, p_state);
		std::vector<CellPtr> cells;

		/*
		 * Create cells as before. this time we use a quiescent volume fraction of x so that
		 * cells
		 */
		for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
		{
			ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
			p_cycle_model->SetCellProliferativeType(TRANSIT);
			p_cycle_model->SetDimension(2);
			p_cycle_model->SetBirthTime(-(double)i - 2.0); // So all out of M phase
			p_cycle_model->SetQuiescentVolumeFraction(0.9);
			p_cycle_model->SetEquilibriumVolume(1.0);

			CellPtr p_cell(new Cell(p_state, p_cycle_model));
			p_cell->InitialiseCellCycleModel();
			cells.push_back(p_cell);
		}

		// Create cell population.
		VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* In order to visualise labelled cells (i.e those that are inhibited from division) you need to use the following command.*/
        cell_population.SetOutputCellMutationStates(true);

        /* We increase the damping constant for healthy cells so the vertices move more slowly */
		cell_population.SetDampingConstantNormal(2*cell_population.GetDampingConstantNormal());

        /* To keep track of the volumes of the cells that are used in the contact inhibition cell-cycle,
         * we use the singleton class {{{CellWiseData}}}. Here, we just initialise it with one variable
         * and associate it with the cell population. This time each cell is associated with a vertex element */
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(1.0, cell_population.GetLocationIndexUsingCell(*cell_iter));
        }

        /*  Then, we define the contact {{{VolumeTrackedOffLatticeSimulation}}} class, that automatically updates the volumes of the cells
         * in {{{CellWiseData}}}. We also set up the output directory, the end time and the output multiple.
         */
		VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
		simulator.SetOutputDirectory("TestVertexContactInhibition");
		simulator.SetSamplingTimestepMultiple(50);
		simulator.SetEndTime(10.0);

		/* Next, we create a force law, {{{NagaiHondaForce}}}, to be applied toe vertices.
		 * We then pass this to the {{{VolumeTrackedOffLatticeSimulation}}} */
		MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* Finally, as in previous cell-based Chaste tutorials, we call {{{Destroy()}}} on the singleton classes. */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellwiseData<2>::Destroy();
	}
    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/TestVertexContactInhibition/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
     * java executable.
     *
     * You will notice that once the healthy cells (yellow) are below a certain size they no longer proliferate and turn dark blue in the visualisation.
     * If you run the simulation for a long time these Cells occur primarily towards the centre of the monolayer.
     *
     * EMPTY LINE
     */

};
#endif /*TESTRUNNIGCONTACTINHIBITIONSIMULATIONSTUTORIAL_HPP_*/
