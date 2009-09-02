/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTCREATINGATISSUESIMULATIONWITHBOUNDARYCONDITIONSTUTORIAL_HPP_
#define TESTCREATINGATISSUESIMULATIONWITHBOUNDARYCONDITIONSTUTORIAL_HPP_

/*
 * = An example showing how to create a tissue simulation with boundary conditions =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how to create a new tissue simulation class in which cells
 * are constrained to lie within a fixed domain.
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


/* The next header defines a base class for tissue simulations, from which
 * the new class will inherit. */
#include "TissueSimulation.hpp"
/* The remaining header files define classes that will be used in the tissue
 * simulation test: {{{HoneycombMeshGenerator}}} defines a helper class for
 * generating a suitable mesh; {{{FixedDurationGenerationBasedCellCycleModelCellsGenerator}}}
 * defines a helper class for generating a vector of cells with fixed cell
 * cycle models; {{{GeneralisedLinearSpringForce}}} defines a force law for
 * describing the mechanical interactions between neighbouring cells in the
 * tissue; and {{{TissueSimulation}}} defines the class that simulates the
 * evolution of a tissue. */
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"

/*
 * EMPTYLINE
 *
 * == Defining the tissue simulation class ==
 *
 * As an example, let us consider a two-dimensional tissue simulation in which
 * all cells are constrained to lie within the domain given in Cartesian coordinates
 * by 0 <= y <= 5. To implement this we define a tissue simulation class,
 * {{{MyTissueSimulation}}}, which inherits from {{{TissueSimulation}}} and overrides
 * the {{{ApplyTissueBoundaryConditions()}}} method.
 */
class MyTissueSimulation : public TissueSimulation<2>
{
/* The first public method is a default constructor, which just calls the base
 * constructor. There are four input argument: a reference to a tissue object,
 * {{{rTissue}}}; a collection of force laws governing tissue mechanics,
 * {{{forceCollection}}}; an optional flag, {{{deleteTissueAndForceCollection}}},
 * telling the simulation whether to delete the tissue and force collection
 * on destruction to free up memory; and another optional flag, {{{initialiseCells}}},
 * telling the simulation whether to initialise cells. */
public:

    MyTissueSimulation(AbstractTissue<2>& rTissue,
                       std::vector<AbstractForce<2>*> forceCollection,
                       bool deleteTissueAndForceCollection=false,
                       bool initialiseCells=true)
        : TissueSimulation<2>(rTissue, forceCollection, deleteTissueAndForceCollection, initialiseCells)
    {
    }

    /* The second public method overrides {{{ApplyTissueBoundaryConditions()}}}.
     * This method is called during the {{{Solve()}}} method at the end of each
     * timestep, just after the position of each node in the tissue has been
     * updated according to its equation of motion. We iterate over all nodes
     * associated with real cells and update their positions according to any
     * boundary conditions. In our case, this means that if any node moves
     * This method iterates over all cells in the tissue, and calls {{{Kill()}}} on
     * any cell whose centre has y coordinate less than 0 or greater than 5. */
    void ApplyTissueBoundaryConditions(const std::vector<c_vector<double,2> >& rOldLocations)
    {
        for (AbstractTissue<2>::Iterator cell_iter = mrTissue.Begin();
             cell_iter != mrTissue.End();
             ++cell_iter)
        {
            c_vector<double, 2> cell_location = mrTissue.GetLocationOfCellCentre(&(*cell_iter));

            unsigned node_index = mrTissue.GetLocationIndexUsingCell(&(*cell_iter));
            Node<2>* p_node = mrTissue.GetNode(node_index);

            if (cell_location[1] > 5.0)
            {
                p_node->rGetModifiableLocation()[1] = 5.0;
            }
            else if (cell_location[1] < 0.0)
            {
                p_node->rGetModifiableLocation()[1] = 0.0;
            }

            assert(p_node->rGetLocation()[1] <= 5.0);
            assert(p_node->rGetLocation()[1] >= 0.0);
        }
    }
};

/* You only need to include the next block of code if you want to be able to
 * archive (save or load) the tissue simulation. We start by including a
 * serialization header, then define {{{save_construct_data}}} and
 * {{{load_construct_data}}} methods, which archive the tissue simulation
 * constructor input argument(s) (in this case, a tissue and a collection of
 * force laws). */
#include "TemplatedExport.hpp"
CHASTE_CLASS_EXPORT(MyTissueSimulation)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const MyTissueSimulation * t, const BOOST_PFTO unsigned int file_version)
        {
            // Save data required to construct instance
            const AbstractTissue<2> * p_tissue = &(t->rGetTissue());
            ar & p_tissue;
            const std::vector<AbstractForce<2>*> force_collection = t->rGetForceCollection();
            ar & force_collection;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, MyTissueSimulation * t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            AbstractTissue<2>* p_tissue;
            ar >> p_tissue;
            std::vector<AbstractForce<2>*> force_collection;
            ar >> force_collection;

            // Invoke inplace constructor to initialise instance
            ::new(t)MyTissueSimulation(*p_tissue, force_collection, true, false);
        }
    }
}


/*
 * This completes the code for {{{MyTissueSimulation}}}. Note that usually this code
 * would be separated out into a separate declaration in a .hpp file and definition
 * in a .cpp file.
 *
 * EMPTYLINE
 *
 * === The Tests ===
 *
 * EMPTYLINE
 *
 * We now define the test class, which inherits from {{{CxxTest::TestSuite}}}.
 */
class TestCreatingATissueSimulationWithBoundaryConditionsTutorial : public CxxTest::TestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Testing the tissue simulation ==
     *
     * EMPTYLINE
     *
     * We now test that our new tissue simulation is implemented correctly.
     */
    void TestMyTissueSimulation() throw(Exception)
    {
        /* The first thing to do, as before, is to set up the start time and
         * reset the parameters. */
        SimulationTime::Instance()->SetStartTime(0.0);
        TissueConfig::Instance()->Reset();

        /* We use the honeycomb mesh generator to create a honeycomb mesh. */
        HoneycombMeshGenerator generator(5, 5, 0, false);
        /* Get the mesh using the {{{GetMesh()}}} method. */
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{TissueCell}}}s.
         * To do this, we can use a static method on the
         * {{{FixedDurationGenerationBasedCellCycleModelCellsGenerator}}} helper class.
         * The {{{<2>}}} below denotes the dimension. We create an empty vector of cells
         * and pass this into the method along with the mesh. The {{{cells}}} vector is
         * populated once the method {{{GenerateBasic}}} is called. */
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        /* Now that we have defined the mesh and cells, we can define the tissue. The
         * constructor takes in the mesh and the cells vector. */
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        /* We must now create one or more force laws, which determine the mechanics of
         * the tissue. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring, and so use
         * a {{{GeneralisedLinearSpringForce}}} object. We pass a pointer to this force
         * into a vector. Note that we have called the method {{{UseCutoffPoint}}} on the
         * {{{GeneralisedLinearSpringForce}}} before passing it into the collection of force
         * laws - this modifies the force law so that two neighbouring cells do not impose
         * a force on each other if they are located more than 3 units (=3 cell widths)
         * away from each other. This modification is necessary when no ghost nodes are used,
         * for example to avoid artificially large forces between cells that lie close together
         * on the tissue boundary.
         */
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(3);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        /* We pass in the tissue and the mechanics system into a {{{TissueSimulation}}}. */
        MyTissueSimulation simulator(tissue, force_collection);

        /* We set the output directory and end time. */
        simulator.SetOutputDirectory("TestMyTissueSimulation");
        simulator.SetEndTime(30.0);

        /* Test that the Solve() method does not throw any exceptions: */
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        /* Test that the boundary conditions have been implemented properly: */
        for (AbstractTissue<2>::Iterator cell_iter = simulator.rGetTissue().Begin();
             cell_iter != simulator.rGetTissue().End();
             ++cell_iter)
        {
            unsigned node_index = simulator.rGetTissue().GetLocationIndexUsingCell(&(*cell_iter));
            Node<2>* p_node = simulator.rGetTissue().GetNode(node_index);

            TS_ASSERT_LESS_THAN_EQUALS(p_node->rGetModifiableLocation()[1], 5.0);
            TS_ASSERT_LESS_THAN_EQUALS(0.0, p_node->rGetModifiableLocation()[1]);
        }

        /* We conclude the test by calling {{{Destroy()}}} on any singleton classes. */
        SimulationTime::Destroy();
    }
};

#endif /*TESTCREATINGATISSUESIMULATIONWITHBOUNDARYCONDITIONSTUTORIAL_HPP_*/
