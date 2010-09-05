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
#ifndef TESTCREATINGANDUSINGANEWCELLKILLERTUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWCELLKILLERTUTORIAL_HPP_

/*
 * = An example showing how to create a new cell killer and use it in a cell-based simulation =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how to create a new cell killer class and how this
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
 * if you want to be able to archive (save or load) the new cell killer object
 * in a cell-based simulation (in this case, these headers must be included before
 * any other serialisation headers). */
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

/* The next header defines a base class for cell killers, from which the new
 * cell killer class will inherit. */
#include "AbstractCellKiller.hpp"
/* The remaining header files define classes that will be used in the cell population
 * simulation test: {{{HoneycombMeshGenerator}}} defines a helper class for
 * generating a suitable mesh; {{{CellsGenerator}}}
 * defines a helper class for generating a vector of cells and
 * {{{FixedDurationGenerationBasedCellCycleModel}}} makes them have fixed cell
 * cycle models; {{{GeneralisedLinearSpringForce}}} defines a force law for
 * describing the mechanical interactions between neighbouring cells in the
 * cell population; and {{{CellBasedSimulation}}} defines the class that simulates the
 * evolution of a cell population. */
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellBasedSimulation.hpp"
#include "CellsGenerator.hpp"
/*
 * EMPTYLINE
 *
 * == Defining the cell killer class ==
 *
 * As an example, let us consider a cell killer which labels any cells in a
 * two-dimensional cell population which lie outside the elliptical domain given in
 * Cartesian coordinates by the equation (x/20)^2^ + (y/10)^2^ < 1. To
 * implement this we define a new cell killer class, {{{MyCellKiller}}},
 * which inherits from {{{AbstractCellKiller}}} and overrides the
 * {{{TestAndLabelCellsForApoptosisOrDeath()}}} method.
 */
class MyCellKiller : public AbstractCellKiller<2>
{
private:

    /* You only need to include the next block of code if you want to be able
     * to archive (save or load) the cell killer object in a cell-based simulation.
     * The code consists of a serialize method, which in this case just archives
     * the cell killer using the serialization code defined in the base class
     * {{{AbstractCellKiller}}}. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
    }

/* The first public method is a default constructor, which just calls the base
 * constructor. */
public:

    MyCellKiller(AbstractCellPopulation<2>* pCellPopulation)
        : AbstractCellKiller<2>(pCellPopulation)
    {}

    /* The second public method overrides {{{TestAndLabelCellsForApoptosisOrDeath()}}}.
     * This method iterates over all cells in the cell_population, and calls {{{Kill()}}} on
     * any cell whose centre is located outside the ellipse (x/20)^2^ + (y/10)^2^ < 1. */
    void TestAndLabelCellsForApoptosisOrDeath()
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
            cell_iter != this->mpCellPopulation->End();
            ++cell_iter)
        {
            c_vector<double, 2> location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

            if ( pow(location[0]/20, 2) + pow(location[1]/10, 2) > 1.0 )
            {
                cell_iter->Kill();
            }
        }
    }
};

/* You only need to include the next block of code if you want to be able to
 * archive (save or load) the cell killer object in a cell-based simulation. We
 * start by including a serialization header, then define {{{save_construct_data}}}
 * and {{{load_construct_data}}} methods, which archive the cell killer
 * constructor input argument(s) (in this case, a {{{CellPopulation}}}). */

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MyCellKiller)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const MyCellKiller * t, const BOOST_PFTO unsigned int file_version)
        {
            // Save data required to construct instance
            const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, MyCellKiller * t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            AbstractCellPopulation<2>* p_cell_population;
            ar >> p_cell_population;

            // Invoke inplace constructor to initialise instance
            ::new(t)MyCellKiller(p_cell_population);
        }
    }
}


/*
 * This completes the code for {{{MyCellKiller}}}. Note that usually this code
 * would be separated out into a separate declaration in a .hpp file and definition
 * in a .cpp file.  In this case, serialization has to be handled slightly more
 * carefully; see BoostSerialization for details.
 *
 * EMPTYLINE
 *
 * === The Tests ===
 *
 * EMPTYLINE
 *
 * We now define the test class, which inherits from {{{CxxTest::TestSuite}}}.
 */
class TestCreatingAndUsingANewCellKillerTutorial : public CxxTest::TestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Testing the cell killer ==
     *
     * EMPTYLINE
     *
     * We begin by testing that our new cell cycle model is implemented correctly.
     */
    void TestMyCellKiller() throw(Exception)
    {
        /* The first thing to do is to set up the start time and reset the model
         * parameters. */
        SimulationTime::Instance()->SetStartTime(0.0);
        CellBasedConfig::Instance()->Reset();

        /* We use the honeycomb mesh generator to create a honeycomb mesh. */
        HoneycombMeshGenerator generator(20, 20, 0, false);
        /* Get the mesh using the {{{GetMesh()}}} method. */
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we can use a static method on the {{{CellsGenerator}}} helper class.
         * The {{{<FixedDurationGenerationBasedCellCycleModel, 2>}}} defines the
         * cell-cycle model and that it is in 2d. We create an empty vector of cells
         * and pass this into the method along with the mesh. The {{{cells}}} vector is
         * populated once the method {{{GenerateBasic}}} is called. */
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        /* Now that we have defined the mesh and cells, we can define the cell population. The
         * constructor takes in the mesh and the cells vector. */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We now use the cell population to construct a cell killer object. */
        MyCellKiller my_cell_killer(&cell_population);

        /* To test that we have implemented the cell killer correctly, we call the
         * overridden method {{{TestAndLabelCellsForApoptosisOrDeath}}}... */
        my_cell_killer.TestAndLabelCellsForApoptosisOrDeath();

        /* ... and check that any cell whose centre is located outside the ellipse
         * (x/20)^2^ + (y/10)^2^ < 1 has indeed been labelled as dead. */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

            if ( pow(x/20, 2) + pow(y/10, 2) > 1.0 )
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
            }
        }

        /* As an extra test, we now remove any dead cells and check that all
         * remaining cells are indeed located within the ellipse. */
        cell_population.RemoveDeadCells();

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

            TS_ASSERT_LESS_THAN_EQUALS(pow(x/20, 2) + pow(y/10, 2) > 1.0, 1.0);
        }

        /* The last chunk of code provides an archiving test for the cell killer.
         * We create an output archive, save the existing cell killer object via
         * a pointer, then create an input archive and load the cell killer. If
         * the cell killer had any member variables, then we would test that these
         * were correctly initialised when the cell killer is loaded. */
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "my_cell_killer.arch";

        {
            // Create an output archive
            MyCellKiller my_cell_killer(NULL);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            MyCellKiller* const p_cell_killer = &my_cell_killer;
            output_arch << p_cell_killer;
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            MyCellKiller* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            delete p_cell_killer;
        }

        /* We conclude the test by calling {{{Destroy()}}} on any singleton classes. */
        SimulationTime::Destroy();
    }

    /*
     * EMPTYLINE
     *
     * == Using the cell killer in a cell-based simulation ==
     *
     * EMPTYLINE
     *
     * We now provide a test demonstrating how {{{MyCellKiller}}} can be used
     * in a cell-based simulation.
     */
    void TestCellBasedSimulationWithMyCellKiller() throw(Exception)
    {
        /* The first thing to do, as before, is to set up the start time and
         * reset the parameters. */
        SimulationTime::Instance()->SetStartTime(0.0);
        CellBasedConfig::Instance()->Reset();

        /* We use the honeycomb mesh generator to create a honeycomb mesh. */
        HoneycombMeshGenerator generator(20, 20, 0, false);
        /* Get the mesh using the {{{GetMesh()}}} method. */
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we can use a static method on the {{{CellsGenerator}}} helper class.
         * The {{{<FixedDurationGenerationBasedCellCycleModel, 2>}}} defines the
         * cell-cycle model and that it is in 2d. We create an empty vector of cells
         * and pass this into the method along with the mesh. The {{{cells}}} vector is
         * populated once the method {{{GenerateBasic}}} is called. */
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        /* Now that we have defined the mesh and cells, we can define the cell population. The
         * constructor takes in the mesh and the cells vector. */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We now use the cell population to construct a cell killer object. */
        MyCellKiller my_cell_killer(&cell_population);

        /* We must now create one or more force laws, which determine the mechanics of
         * the cell population. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring, and so use
         * a {{{GeneralisedLinearSpringForce}}} object. We pass a pointer to this force
         * into a vector. Note that we have called the method {{{UseCutoffPoint}}} on the
         * {{{GeneralisedLinearSpringForce}}} before passing it into the collection of force
         * laws - this modifies the force law so that two neighbouring cells do not impose
         * a force on each other if they are located more than 3 units (=3 cell widths)
         * away from each other. This modification is necessary when no ghost nodes are used,
         * for example to avoid artificially large forces between cells that lie close together
         * on the cell population boundary.
         */
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(3);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        /* We pass in the cell population and the mechanics system into a {{{CellBasedSimulation}}}. */
        CellBasedSimulation<2> simulator(cell_population, force_collection);

        /* We set the output directory and end time. */
        simulator.SetOutputDirectory("TestCellBasedSimulationWithMyCellKiller");
        simulator.SetEndTime(10.0);

        /* We now pass the cell killer into the cell-based simulation. */
        MyCellKiller* p_killer = new MyCellKiller(&cell_population);
        simulator.AddCellKiller(p_killer);

        /* Test that the Solve() method does not throw any exceptions. */
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        /* We conclude the test by calling {{{Destroy()}}} on any singleton classes. */
        SimulationTime::Destroy();
    }
};

#endif /*TESTCREATINGANDUSINGANEWCELLKILLERTUTORIAL_HPP_*/
