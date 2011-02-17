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
#ifndef TESTCREATINGANDUSINGANEWCELLCYCLEMODELTUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWCELLCYCLEMODELTUTORIAL_HPP_

/*
 * = An example showing how to create a new cell cycle model and use it in a cell-based simulation =
 *
 * == Introduction ==
 *
 * In this tutorial we show how to create a new cell cycle model class and how this
 * can be used in a cell-based simulation.
 *
 * == 1. Including header files ==
 *
 * The first thing to do is include the following header, which allows us
 * to use certain methods in our test (this header file should be included
 * in any Chaste test):
 */
#include <cxxtest/TestSuite.h>

/* The next two headers are used in archiving, and only need to be included
 * if you want to be able to archive (save or load) the new cell killer object
 * in a cell-based simulation (in this case, these headers must be included before
 * any other serialisation headers). */
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

/* The next header defines a base class for simple generation-based cell
 * cycle models.
 *
 * A cell cycle model is 'simple' if the duration of each phase of the cell
 * cycle is determined when the cell cycle model is created, rather than
 * evaluated 'on the fly' (e.g. by solving a system of ordinary differential
 * equations for the concentrations of key cell cycle proteins), and may
 * depend on the cell type.
 *
 * A simple cell cycle model is generation-based if it keeps track of the
 * generation of the corresponding cell, and sets the cell type according
 * to this.
 *
 * Our new cell cycle model will inherit from this abstract class. */
#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"

/* The remaining header files define classes that will be used in the cell population
 * simulation test: {{{CheckReadyToDivideAndPhaseIsUpdated}}} defines a helper
 * class for testing a cell cycle model; {{{HoneycombMeshGenerator}}} defines
 * a helper class for generating a suitable mesh; {{{WildTypeCellMutationState}}}
 * defines a wild-type or 'healthy' cell mutation state; {{{GeneralisedLinearSpringForce}}}
 * defines a force law for describing the mechanical interactions between neighbouring
 * cells in the cell population; and {{{CellBasedSimulation}}} defines the class that
 * simulates the evolution of the cell population. */
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellBasedSimulation.hpp"

/*
 * == Defining the cell cycle model class ==
 *
 * As an example, let us consider a cell cycle model in which the durations
 * of S, G2 and M phases are fixed, but the duration of G1 phase is an exponential
 * random variable with rate parameter lambda.
 *
 * The rate parameter is a constant, dependent on cell type, whose value is
 * chosen such that the mean of the distribution, 1/lambda, equals the mean
 * G1 duration as defined in the {{{AbstractCellCycleModel}}} class.
 *
 * To implement this model we define a new cell cycle model, {{{MyCellCycleModel}}},
 * which inherits from {{{AbstractSimpleGenerationBasedCellCycleModel}}} and
 * overrides the {{{SetG1Duration()}}} method.
 * 
 * Note that usually this code would be separated out into a separate declaration in
 * a .hpp file and definition in a .cpp file.
 */
class MyCellCycleModel : public AbstractSimpleGenerationBasedCellCycleModel
{
private:

    /* You only need to include the next block of code if you want to be able
     * to archive (save or load) the cell cycle model object in a cell-based simulation.
     * The code consists of a serialize method, in which we first archive the cell
     * cycle model using the serialization code defined in the base class
     * {{{AbstractSimpleGenerationBasedCellCycleModel}}}. We then archive an instance
     * of the {{{RandomNumberGenerator}}} singleton class, which is used in the
     * {{{SetG1Duration()}}} method. Note that serialization of singleton objects
     * must be done with care. Before the object is serialized via a pointer, it must
     * be serialized directly, or an assertion will trip when a second instance of the
     * class is created on de-serialization. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleGenerationBasedCellCycleModel>(*this);
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
    }

    /* We override the {{{SetG1Duration()}}} method as follows. */
    void SetG1Duration()
    {
        /* As we will access the cell type of the cell associated with this cell
         * cycle model, we should assert that this cell exists. */
        assert(mpCell != NULL);

        /* We now set the G1 duration based on cell type.
         *
         * For stem and transit cells, we use the {{{RandomNumberGenerator}}}
         * singleton class to generate a random number U drawn from U![0,1], and
         * transform this into a random number T drawn from Exp(lambda) using
         * the transformation T = -log(U)/lambda.
         *
         * For differentiated cells, which do not progress through the
         * cell cycle, we set the G1 duration to {{{DBL_MAX}}}. */
        double uniform_random_number = RandomNumberGenerator::Instance()->ranf();

        switch (mCellProliferativeType)
        {
            case STEM:
                mG1Duration = -log(uniform_random_number)*GetStemCellG1Duration();
                break;
            case TRANSIT:
                mG1Duration = -log(uniform_random_number)*GetTransitCellG1Duration();
                break;
            case DIFFERENTIATED:
                mG1Duration = DBL_MAX;
                break;
            default:
                NEVER_REACHED;
        }
    }

/* The first public method is a default constructor, which just calls the base
 * constructor. */
public:

    MyCellCycleModel()
    {}

    /* The second public method overrides {{{CreateCellCycleModel()}}}. This is a
     * builder method to create new copies of the cell cycle model. */
    AbstractCellCycleModel* CreateCellCycleModel()
    {
        // Create a new cell cycle model
        MyCellCycleModel* p_model = new MyCellCycleModel();
        
        // Set the values of the new cell cycle model's member variables
        p_model->SetGeneration(mGeneration);
        p_model->SetMaxTransitGenerations(mMaxTransitGenerations);

        return p_model;
    }
};

/* You need to include the next block of code if you want to be able to archive (save or load)
 * the cell cycle model object in a cell-based simulation.  It is also required for writing out
 * the parameters file describing the settings for a simulation - it provides the unique
 * identifier for our new cell cycle model.  Thus every cell cycle model class must provide this,
 * or you'll get errors when running simulations. */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MyCellCycleModel)

/* Since we're defining the new cell cycle model within the test file, we need to include the
 * following stanza as well, to make the code work with newer versions of the Boost libraries.
 * Normally the above export declaration would occur in the cell cycle model's .hpp file, and
 * the following lines would appear in the .cpp file.  See ChasteGuides/BoostSerialization for
 * more information.
 */
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyCellCycleModel)

/*
 * This completes the code for {{{MyCellCycleModel}}}. Note that usually this code would
 * be separated out into a separate declaration in a .hpp file and definition in a .cpp file.
 *
 * === The Tests ===
 *
 * We now define the test class, which inherits from {{{CxxTest::TestSuite}}}.
 */
class TestCreatingAndUsingANewCellCycleModelTutorial : public CxxTest::TestSuite
{
public:

    /*
     * == Testing the cell cycle model ==
     *
     * We begin by testing that our new cell cycle model is implemented correctly.
     */
    void TestMyCellCycleModel() throw(Exception)
    {
        /* We must first set the start time. In addition, it is advisable to reset
         * the values of all model parameters. Recall that {{{SimulationTime}}} is a 
         * ''singleton'' class; this means one and only
         * one of each of this object is instantiated at any time, and that single
         * object is accessible from anywhere in the code. As a result, we do not need
         * to keep passing round the current time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Test that we can construct a {{{MyCellCycleModel}}} object: */
        TS_ASSERT_THROWS_NOTHING(MyCellCycleModel cell_model3);

        /* Now construct and initialise a large number of {{{MyCellCycleModel}}}s and
         * associated cells: */
        unsigned num_cells = 1e5;
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<num_cells; i++)
        {
            MyCellCycleModel* p_cell_cycle_model = new MyCellCycleModel;
            p_cell_cycle_model->SetCellProliferativeType(STEM);
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        /* Find the mean G1 duration and test that it is within some tolerance of
         * the expected value: */
        double expected_mean_g1_duration = cells[0]->GetCellCycleModel()->GetStemCellG1Duration();
        double sample_mean_g1_duration = 0.0;

        for (unsigned i=0; i<num_cells; i++)
        {
            sample_mean_g1_duration += cells[i]->GetCellCycleModel()->GetG1Duration()/ (double) num_cells;
        }

        TS_ASSERT_DELTA(sample_mean_g1_duration, expected_mean_g1_duration, 0.1);

        /* Now construct another {{{MyCellCycleModel}}} and associated cell. */
        MyCellCycleModel* p_my_model = new MyCellCycleModel;
        p_my_model->SetCellProliferativeType(TRANSIT);
        CellPtr p_my_cell(new Cell(p_state, p_my_model));
        p_my_cell->InitialiseCellCycleModel();

        /* Use the helper method {{{CheckReadyToDivideAndPhaseIsUpdated()}}} to
         * test that this cell progresses correctly through the cell cycle. */
        unsigned num_steps = 100;
        double mean_cell_cycle_time = cells[0]->GetCellCycleModel()->GetStemCellG1Duration()
                                        + cells[0]->GetCellCycleModel()->GetSG2MDuration();

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(mean_cell_cycle_time, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();

            /* The numbers for the G1 duration below is taken from the first
             * random number generated: */
            CheckReadyToDivideAndPhaseIsUpdated(p_my_model, 1.18892);
        }

        /* Lastly, we briefly test that archiving of {{{MyCellCycleModel}}} has
         * been implemented correctly. Create an {{{OutputFileHandler}}} and use
         * this to define a filename for the archive.
         */
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "my_cell_cycle_model.arch";

        /* Create an output archive. */
        {
            /* Destroy the current instance of {{{SimulationTime}}} and create another instance.
             * Set the start time, end time and number of time steps. */
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 4);

            /* Create a cell with associated cell cycle model. */
            MyCellCycleModel* p_model = new MyCellCycleModel;
            p_model->SetCellProliferativeType(TRANSIT);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->InitialiseCellCycleModel();

            /* Move forward two time steps. */
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();

            /* Set the birth time of the cell and update the cell cycle phase. */
            p_model->SetBirthTime(-1.0);
            p_model->ReadyToDivide();

            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), S_PHASE);

            /* Now archive the cell cycle model through its cell. */
            CellPtr const p_const_cell = p_cell;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_const_cell;
        }

        /* Now create an input archive. Begin by again destroying the current
         * instance of {{{SimulationTime}}} and creating another instance. Set
         * the start time, end time and number of time steps.
         */
        {
            SimulationTime::Destroy();
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            /* Create a pointer to a cell. */
            CellPtr p_cell;

            /* Create an input archive and restore the cell from the archive. */
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell;

            /* Test that the private data has been restored correctly. */
            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();

            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(), 2.5, 1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), S_PHASE);
        }

        /* {{{SimulationTime::Destroy()}}} '''must''' be called at the end of the test.
         * If not, when {{{SimulationTime::Instance()->SetStartTime(0.0);}}} is called
         * at the beginning of the next test in this file, an assertion will be triggered. */
        SimulationTime::Destroy();
        /* Also call {{{Destroy()}}} on the {{{RandomNumberGenerator}}} singleton class. */
        RandomNumberGenerator::Destroy();
    }

    /*
     * == Using the cell cycle model in a cell-based simulation ==
     *
     * We conclude with a brief test demonstrating how {{{MyCellCycleModel}}} can be used
     * in a cell-based simulation.
     */
    void TestCellBasedSimulationWithMyCellCycleModel() throw(Exception)
    {
        /* The first thing to do, as before, is to set up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* We use the honeycomb mesh generator to create a honeycomb mesh covering a
         * circular domain of given radius.
         */
        HoneycombMeshGenerator generator(10, 10, 0);
        /* Get the mesh using the {{{GetCircularMesh()}}} method. */
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(5);

        /* Next, we create some cells. First, define the cells vector. */
        std::vector<CellPtr> cells;
        /* Then we loop over the nodes. */
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            /* For each node we create a cell with our cell cycle model. */
            MyCellCycleModel* p_model = new MyCellCycleModel();
            p_model->SetCellProliferativeType(STEM);
            CellPtr p_cell(new Cell(p_state, p_model));

            /* Now, we define a random birth time, chosen from [-T,0], where
             * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
             * of a stem cell, and t,,2,, is the basic S+G,,2,,+M phases duration.
             */
            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                    (p_model->GetStemCellG1Duration()
                                        + p_model->GetSG2MDuration());
            /* We then set the birth time and push the cell back into the vector of cells. */
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Now that we have defined the mesh and cells, we can define the cell population. The constructor
         * takes in the mesh and the cells vector. */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /*
         * We pass in the cell population into a {{{CellBasedSimulation}}}.
         */
        CellBasedSimulation<2> simulator(cell_population);

        /* We must now create one or more force laws, which determine the mechanics of
         * the cell population. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring, and so use
         * a {{{GeneralisedLinearSpringForce}}} object. We pass a pointer to this force
         * into a vector. Note that we have called the method {{{SetCutOffLength}}} on the
         * {{{GeneralisedLinearSpringForce}}} before passing it into the collection of force
         * laws - this modifies the force law so that two neighbouring cells do not impose
         * a force on each other if they are located more than 3 units (=3 cell widths)
         * away from each other. This modification is necessary when no ghost nodes are used,
         * for example to avoid artificially large forces between cells that lie close together
         * on the cell population boundary.
         */

        /* We create a force law and pass it to the {{{CellBasedSimulation}}}. */
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(3);
        simulator.AddForce(&linear_force);

        /* We set the output directory and end time. */
        simulator.SetOutputDirectory("TestCellBasedSimulationWithMyCellCycleModel");
        simulator.SetEndTime(10.0);

        /* Test that the Solve() method does not throw any exceptions. */
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        /* Finally, call {{{Destroy()}}} on the singleton classes. */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTCREATINGANDUSINGANEWCELLCYCLEMODELTUTORIAL_HPP_*/
