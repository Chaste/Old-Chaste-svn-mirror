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
#ifndef TESTCREATINGANDUSINGANEWCELLMUTATIONSTATETUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWCELLMUTATIONSTATETUTORIAL_HPP_

/*
 * = An example showing how to create a new cell mutation state and use it in a cell-based simulation =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how to create a new cell mutation state class and how this
 * can be used in a cell-based simulation.
 *
 * EMPTYLINE
 *
 * == 1. Including header files ==
 *
 * EMPTYLINE
 *
 * The first thing to do is include the following header, which allows us
 * to use certain methods in our test (this header file should be included
 * in any Chaste test):
 */
#include <cxxtest/TestSuite.h>

/* The next two headers are used in archiving, and only need to be included
 * if we intend to archive (save or load) a cell-based simulation in this test
 * suite. In this case, these headers must be included before any other
 * serialisation headers. To use archiving uncomment these lines.
 * Note, on machines running boost 1.40 you need to split the code between a
 * cpp and hpp file for archiving to work, therefore it is commented here. */
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>

/* The next header defines a base class for cell mutation states. Our new
 * cell mutation state will inherit from this abstract class. */
#include "AbstractCellMutationState.hpp"
/* The remaining header files define classes that will be used in the cell population
 * simulation test: {{{HoneycombMeshGenerator}}} defines a helper class for
 * generating a suitable mesh; {{{WildTypeCellMutationState}}} defines a
 * wild-type or 'healthy' cell mutation state; {{{FixedDurationGenerationBasedCellCycleModel}}}
 * defines a simple cell-cycle model class, in which cells undergo a fixed number
 * of divisions before becoming senescent; {{{GeneralisedLinearSpringForce}}}
 * defines a force law for describing the mechanical interactions between neighbouring
 * cells in the cell population; and {{{CellBasedSimulation}}} defines the class that
 * simulates the evolution of the cell population. */
#include "HoneycombMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellBasedSimulation.hpp"
#include "CellsGenerator.hpp"

/*
 * EMPTYLINE
 *
 * == Defining the cell mutation state class ==
 *
 * As an example, let us consider a cell mutation state representing the p53
 * 172R-H gain-of-function mutant, which is equivalent to the common 175R-H
 * human breast cancer mutant; for further details on this mutant, see for
 * example Murphy et al, FASEB J. 14:2291-2302 (2000).
 *
 * Wild-type p53 has been referred to as the "guardian of the genome",
 * responding to DNA damage or checkpoint failure by either arresting cell
 * cycle progression to facilitate DNA repair or initiating an apoptotic
 * pathway to remove damaged cells. Approximately 40% of human breast cancers
 * contain alterations in p53.
 *
 * As we can see, apart from a serialize() method and a constructor, this class
 * does not contain any member variables or methods. This is because generally
 * a cell's mutation state is used, much like a flag, by other classes when
 * determining a cell's behaviour (whether a cell should undergo
 * apoptosis following prolonged stress, for example, or alter its proliferative
 * behaviour).
 */
class P53GainOfFunctionCellMutationState : public AbstractCellMutationState
{
private:
    /* The next block of code allows us to archive (save or load) the cell mutation
     * state object in a cell-based simulation. The code consists of a serialize() method,
     * in which we archive the cell mutation state using the serialization code defined in the base class
     * {{{AbstractCellMutationState}}}. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    /* The only public method is a default constructor, which just calls the base
     * constructor with a single unsigned parameter. This sets the value of the
     * base class member variable {{{mColour}}}, which can be used by visualization tools
     * to paint cells with this mutation state a distinct colour if required. */
    P53GainOfFunctionCellMutationState()
        : AbstractCellMutationState(5)
    {
    }
};

/* Together with the serialize() method defined within the class above, the next
 * block of code allows us to archive (save or load) the cell mutation state
 * object in a cell-based simulation. To use archiving uncomment these lines.
 * Note, on machines running boost 1.40 you need to split the code between a
 * cpp and hpp file for archiving to work, therefore it is commented here. */
//#include "SerializationExportWrapper.hpp"
//CHASTE_CLASS_EXPORT(P53GainOfFunctionCellMutationState)

/*
 * This completes the code for {{{P53GainOfFunctionCellMutationState}}}. Note that usually this code would
 * be separated out into a separate declaration in a .hpp file and definition in a .cpp file.
 *
 * EMPTYLINE
 *
 * === The Tests ===
 *
 * EMPTYLINE
 *
 * We now define the test class, which inherits from {{{CxxTest::TestSuite}}}.
 */
class TestCreatingAndUsingANewCellMutationStateTutorial : public CxxTest::TestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Testing the cell mutation state ==
     *
     * EMPTYLINE
     *
     * We begin by testing that our new cell mutation state is implemented correctly.
     */
    void TestP53GainOfFunctionCellMutationState() throw(Exception)
    {
        /* We begin by testing that some of the base class methods work correctly.
         * We typically use shared pointers to create and access cell mutation states, as
         * follows. This is because it makes sense for all cells that have the same mutation
         * to share a pointer to the same cell mutation state object (although strictly speaking,
         * they are not required to).*/
        boost::shared_ptr<AbstractCellMutationState> p_state(new P53GainOfFunctionCellMutationState);

        /* Each cell mutation state has a member variable, {{{mCellCount}}}, which
         * stores the number of cells with this mutation state. In fact, {{{mCellCount}}}
         * is defined in the class {{{AbstractCellProperty}}}, from which
         * {{{AbstractCellMutationState}}} inherits, as well as other cell properties
         * such as {{{CellLabel}}}. We can test whether {{{mCellCount}}} is being
         * updated correctly by our cell mutation state, as follows. */
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
        p_state->IncrementCellCount();
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
        p_state->DecrementCellCount();
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
        TS_ASSERT_THROWS_THIS(p_state->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");

        /* We can also test that {{{mColour}}} has been set correctly by our constructor, as follows. */
        TS_ASSERT_EQUALS(p_state->GetColour(), 5u);

        /* We can also test whether our cell mutation state is of a given type, as follows. */
        TS_ASSERT_EQUALS(p_state->IsType<WildTypeCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_state->IsType<P53GainOfFunctionCellMutationState>(), true);

        /* We can also test that archiving is implemented correctly for our cell
         * mutation state, as follows (further details on how to implement and
         * test archiving can be found on the BoostSerialization page).
         * To test archiving uncomment these lines. Note, on machines running
         * boost 1.40 you need to split the code between a cpp and hpp file for
         * archiving to work, therefore it is commented here.*/
//        OutputFileHandler handler("archive", false);
//        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "mutation.arch";
//
//        {
//            P53GainOfFunctionCellMutationState* p_state = new P53GainOfFunctionCellMutationState();
//            p_state->IncrementCellCount();
//
//            TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
//            TS_ASSERT_EQUALS(p_state->GetColour(), 5u);
//
//            std::ofstream ofs(archive_filename.c_str());
//            boost::archive::text_oarchive output_arch(ofs);
//
//            const AbstractCellProperty* const p_const_state = p_state;
//            output_arch << p_const_state;
//
//            delete p_state;
//        }
//
//        {
//            AbstractCellProperty* p_state;
//
//            std::ifstream ifs(archive_filename.c_str());
//            boost::archive::text_iarchive input_arch(ifs);
//
//            input_arch >> p_state;
//
//            TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
//
//            P53GainOfFunctionCellMutationState* p_real_state = dynamic_cast<P53GainOfFunctionCellMutationState*>(p_state);
//            TS_ASSERT(p_real_state != NULL);
//            TS_ASSERT_EQUALS(p_real_state->GetColour(), 5u);
//
//            delete p_state;
//        }
    }

    /*
     * EMPTYLINE
     *
     * == Using the cell mutation state in a cell-based simulation ==
     *
     * EMPTYLINE
     *
     * We conclude with a brief test demonstrating how {{{P53GainOfFunctionCellMutationState}}} can be used
     * in a cell-based simulation.
     */
    void TestCellBasedSimulationWithP53GainOfFunctionCellMutationState() throw(Exception)
    {
        /* We begin by setting up the start time, as follows. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* We use the {{{HoneycombMeshGenerator}}} to create a honeycomb mesh covering a
         * circular domain of given radius, as follows. */
        HoneycombMeshGenerator generator(10, 10, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(5);

        /* We now create a shared pointer to our new cell mutation state, as follows. */
        boost::shared_ptr<AbstractCellMutationState> p_state(new P53GainOfFunctionCellMutationState);

        /* Next, we create some cells, as follows. */
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        /* Now that we have defined the mesh and cells, we can define the cell population. The constructor
         * takes in the mesh and the cells vector. */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /*
         * We pass in the cell population into a {{{CellBasedSimulation}}}.
         */
        CellBasedSimulation<2> simulator(cell_population);

        /* We set the output directory and end time. */
        simulator.SetOutputDirectory("TestCellBasedSimulationWithp_motile");
        simulator.SetEndTime(10.0);

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

        /* Test that the Solve() method does not throw any exceptions. */
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        /* Finally, call {{{Destroy()}}} on the singleton classes. */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTCREATINGANDUSINGANEWCELLMUTATIONSTATETUTORIAL_HPP_*/
