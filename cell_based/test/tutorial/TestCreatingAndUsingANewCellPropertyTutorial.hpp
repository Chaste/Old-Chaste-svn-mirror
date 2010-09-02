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
#ifndef TESTCREATINGANDUSINGANEWCELLPROPERTYTUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWCELLPROPERTYTUTORIAL_HPP_

/*
 * = An example showing how to create a new cell property and use it in a tissue simulation =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how to create a new cell property class and how this
 * can be used in a tissue simulation.
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
 * if you intend to archive (save or load) a tissue simulation in this test
 * suite. In this case, these headers must be included before any other
 * serialisation headers. */
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

/* The next header defines a base class for cell properties. Our new
 * cell property will inherit from this abstract class. */
#include "AbstractCellProperty.hpp"
/* The remaining header files define classes that will be used in the tissue
 * simulation test: {{{HoneycombMeshGenerator}}} defines a helper class for
 * generating a suitable mesh; {{{WildTypeCellMutationState}}} defines a
 * wild-type or 'healthy' cell mutation state; {{{FixedDurationGenerationBasedCellCycleModel}}}
 * defines a simple cell-cycle model class, in which cells undergo a fixed number
 * of divisions before becoming senescent; {{{GeneralisedLinearSpringForce}}}
 * defines a force law for describing the mechanical interactions between neighbouring
 * cells in the tissue; and {{{TissueSimulation}}} defines the class that
 * simulates the evolution of the tissue. */
#include "HoneycombMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "TissueSimulation.hpp"

/*
 * EMPTYLINE
 *
 * == Defining the cell property class ==
 *
 * As an example, let us consider a cell property class that is used to label
 * those cells that are "motile". This cell property could then be used when
 * implementing some form of chemotaxis down an imposed chemoattractant gradient,
 * as occurs for example when macrophages migrate within a tumour towards high
 * concentrations of the vascular endothelial growth factor VEGF; for further
 * details, see for example Owen et al, J. Theor. Biol.
 * 226: 377-391 (2004).
 */
class MotileCellProperty : public AbstractCellProperty
{
private:

    /* We define a member variable {{{mColour}}}, which can be used by visualization tools
     * to paint cells with this mutation state a distinct colour if required. */
    unsigned mColour;

    /* The next block of code allows us to archive (save or load) the cell property object
     * in a tissue simulation. The code consists of a serialize() method, in which we first
     * archive the cell property using the serialization code defined in the base class
     * {{{AbstractCellProperty}}}, then archive the member variable {{{mColour}}}. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

public:

    /* The default constructor allows us to specify a value for the member variable {{{mColour}}},
     * or leave it with a default value. */
    MotileCellProperty(unsigned colour=5)
        : AbstractCellProperty(),
          mColour(colour)
    {
    }

    /* We then define a destructor and a get method for the member variable {{{mColour}}}. */
    ~MotileCellProperty()
    {}

    unsigned GetColour() const
    {
        return mColour;
    }
};

/* Together with the serialize() method defined within the class above, the next
 * block of code allows you to archive (save or load) the cell property object
 * in a tissue simulation. */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MotileCellProperty)

/* This completes the code for {{{MotileCellProperty}}}. Note that usually this code would
 * be separated out into a separate declaration in a .hpp file and definition in a .cpp file.
 * 
 * === The Tests ===
 *
 * EMPTYLINE
 *
 * We now define the test class, which inherits from {{{CxxTest::TestSuite}}}.
 */
class TestCreatingAndUsingANewCellPropertyTutorial : public CxxTest::TestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Testing the cell property ==
     *
     * EMPTYLINE
     *
     * We begin by testing that our new cell property is implemented correctly.
     */
    void TestMotileCellProperty() throw(Exception)
    {
        /* We begin by testing that some of the base class methods work correctly.
         * We typically use shared pointers to create and access a cell property
         * like {{{MotileCellProperty}}}, for which it makes sense for all cells
         * that have the same mutation to share a pointer to the same cell property
         * object (although strictly speaking, they are not required to). Observe that
         * in this case we have provided a value for the member variable {{{mColour}}}
         * in the {{{MotileCellProperty}}} constructor.*/
        boost::shared_ptr<AbstractCellProperty> p_property(new MotileCellProperty(8));

        /* Each cell property has a member variable, {{{mCellCount}}}, which
         * stores the number of cells with this cell property. We can test whether
         * {{{mCellCount}}} is being updated correctly by our cell property, as follows. */
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 0u);
        p_property->IncrementCellCount();
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 1u);
        p_property->DecrementCellCount();
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 0u);
        TS_ASSERT_THROWS_THIS(p_property->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");

        /* We can also test whether our cell property is of a given type, as follows. */
        TS_ASSERT_EQUALS(p_property->IsType<WildTypeCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_property->IsType<MotileCellProperty>(), true);

        /* We can also test that archiving is implemented correctly for our cell
         * property, as follows (further details on how to implement and
         * test archiving can be found on the BoostSerialization page).  */ 
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "property.arch";

        {
            MotileCellProperty* p_property = new MotileCellProperty(7);
            p_property->IncrementCellCount();

            TS_ASSERT_EQUALS(p_property->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_property->GetColour(), 7u);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            const AbstractCellProperty* const p_const_property = p_property;
            output_arch << p_const_property;

            delete p_property;
        }

        {
            AbstractCellProperty* p_property;

            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_property;

            TS_ASSERT_EQUALS(p_property->GetCellCount(), 1u);

            MotileCellProperty* p_real_property = dynamic_cast<MotileCellProperty*>(p_property);
            TS_ASSERT(p_real_property != NULL);
            TS_ASSERT_EQUALS(p_real_property->GetColour(), 7u);

            delete p_property;
        }
    }

    /*
     * EMPTYLINE
     *
     * == Using the cell property in a tissue simulation ==
     *
     * EMPTYLINE
     *
     * We conclude with a brief test demonstrating how {{{MotileCellProperty}}} can be used
     * in a tissue simulation.
     */
    void TestTissueSimulationWithP53GainOfFunctionCellMutationState() throw(Exception)
    {
        /* We begin by setting up the start time, as follows. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* We use the {{{HoneycombMeshGenerator}}} to create a honeycomb mesh covering a
         * circular domain of given radius, as follows. */
        HoneycombMeshGenerator generator(10, 10, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(5);

        /* We now create a shared pointer to our new property, as follows. */ 
        boost::shared_ptr<AbstractCellProperty> p_motile(new MotileCellProperty);

        /* Next, we create some cells, as follows. */
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        std::vector<TissueCellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            /* For each node we create a cell with our cell cycle model and the wild-type cell mutation state.
             * We then add the property {{{MotileCellProperty}}} to a random selection of the cells, as follows. */
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            CellPropertyCollection collection;
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                collection.AddProperty(p_motile);
            }

            TissueCellPtr p_cell(new TissueCell(p_state, p_model, false, collection));

            /* Now, we define a random birth time, chosen from [-T,0], where
             * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
             * of a stem cell, and t,,2,, is the basic S+G,,2,,+M phases duration.
             */
            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                    (TissueConfig::Instance()->GetStemCellG1Duration()
                                        + TissueConfig::Instance()->GetSG2MDuration());

            /* Finally, we set the birth time and push the cell back into the vector of cells. */
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Now that we have defined the mesh and cells, we can define the tissue. The constructor
         * takes in the mesh and the cells vector. */
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        /* We must now create one or more force laws, which determine the mechanics of
         * the tissue. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring, and so use
         * a {{{GeneralisedLinearSpringForce}}} object. Note that we have called the method {{{UseCutoffPoint}}} on the
         * {{{GeneralisedLinearSpringForce}}} before passing it into the collection of force
         * laws - this modifies the force law so that two neighbouring cells do not impose
         * a force on each other if they are located more than 3 units (=3 cell widths)
         * away from each other. This modification is necessary when no ghost nodes are used,
         * for example to avoid artificially large forces between cells that lie close together
         * on the tissue boundary.
         */
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(3);

        /* We then pass a pointer to the force into a vector. */
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        /*
         * We pass in the tissue and the mechanics system into a {{{TissueSimulation}}}.
         */
        TissueSimulation<2> simulator(tissue, force_collection);

        /* We set the output directory and end time. */
        simulator.SetOutputDirectory("TestTissueSimulationWithMotileCellProperty");
        simulator.SetEndTime(10.0);

        /* Test that the Solve() method does not throw any exceptions. */
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        /* Finally, call {{{Destroy()}}} on the singleton classes. */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTCREATINGANDUSINGANEWCELLPROPERTYTUTORIAL_HPP_*/
