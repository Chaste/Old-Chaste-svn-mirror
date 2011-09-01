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

#ifndef TESTCREATINGANDUSINGANEWFORCETUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWFORCETUTORIAL_HPP_

/*
 * = An example showing how to create and use a new force =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In previous cell-based Chaste tutorial, we used existing force classes to define
 * how cells interact mechanically. In this tutorial we show
 * how to create a new force class, and how this can be used in a cell-based
 * simulation.
 *
 * EMPTYLINE
 *
 * == 1. Including header files ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
/* The next header defines a base class for forces, from which the new class will inherit. */
#include "AbstractForce.hpp"
/* The remaining header files define classes that will be used in the cell population
 * simulation test. We have encountered each of these header files in previous cell-based
 * Chaste tutorials. */
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"

/*
 * EMPTYLINE
 *
 * == Defining the force class ==
 *
 * As an example, let us consider a force for a two-dimensional cell-based
 * simulation, that mimics gravity. To implement this we define a force
 * boundary condition class, {{{MyBoundaryCondition}}}, which inherits from
 * {{{AbstractForce}}} and overrides the methods {{{AddForceContribution()}}} and
 * {{{OutputForceParameters()}}}.
 */
class MyForce : public AbstractForce<2>
{
private:

    /* This force class includes a member variable, {{{mStrength}}}, which
     * defines the strength of the force. This member variables will be set
     * in the constructor.
     */
    double mStrength;

    /* We only need to include the next block of code if we wish to be able
     * to archive (save or load) the force model object in a cell-based simulation.
     * The code consists of a serialize method, in which we first archive the force
     * using the serialization code defined in the base class {{{AbstractForce}}},
     * then archive the member variable. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<2> >(*this);
        archive & mStrength;
    }

public:
    /* The first public method is a default constructor, which calls the base
     * constructor. There is a single input argument, which defines the strength
     * of the force. We provide a default value of 1.0 for this argument. Inside
     * the method, we add an assertion to make sure that the strength is strictly
     * positive.
     */
    MyForce(double strength=1.0)
        : AbstractForce<2>(),
          mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    /* The second public method overrides {{{AddForceContribution()}}}.
     * This method takes in two arguments: a reference to a vector of
     * total forces on nodes in a cell population, which is update to by the
     * force object; and a reference to the cell population itself.
     */
    void AddForceContribution(std::vector<c_vector<double, 2> >& rForces,
                              AbstractCellPopulation<2>& rCellPopulation)
    {
        /* Inside the method, we loop over nodes, and add a constant vector to
         * each node, in the negative ''y''-direction and of magnitude {{{mStrength}}}.
         */
        c_vector<double, 2> force = zero_vector<double>(2);
        force(1) = -mStrength;

        for (unsigned node_index=0; node_index<rForces.size(); node_index++)
        {
            rForces[node_index] += force;
        }
    }

    /* We also add a get method for {{{mStrength}}}, to allow for testing. */
    double GetStrength()
    {
        return mStrength;
    }

    /* Just as we encountered in the cell killer tutorial, here we must override
     * a method that outputs any member variables to a specified results file {{{rParamsFile}}}.
     * In our case, we output the member variable {{{mStrength}, then call the method on the base class.
     */
    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength> \n";
        AbstractForce<2>::OutputForceParameters(rParamsFile);
    }
};

/* As mentioned in previous cell-based Chaste tutorials, we need to include the next block
 * of code to be able to archive the force object in a cell-based
 * simulation, and to obtain a unique identifier for our new force for writing
 * results to file.
 */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MyForce)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyForce)

/*
 * This completes the code for {{{MyForce}}}. Note that usually this code
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
class TestCreatingAndUsingANewForceTutorial : public CxxTest::TestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Testing the force ==
     *
     * EMPTYLINE
     *
     * We now test that our new force is implemented correctly.
     */
    void TestMyForce() throw(Exception)
    {
        /* The first thing to do, as before, is to set up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* We now create a {{{MeshBasedCellPopulation}}} using the helper
         * classes {{{HoneycombMeshGenerator}}} and {{{CellsGenerator}}},
         * as in previous cell-based Chaste tutorials.
         */
        HoneycombMeshGenerator generator(7, 7, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* The next step is to initialise a vector of node forces. */
        std::vector<c_vector<double,2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }

        /* We now create a force object of strength 5.0.
         */
        MyForce force(5.0);

        /* We test that the force calculation is correct. */
        force.AddForceContribution(node_forces, cell_population);

        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(node_forces[node_index][0], 0.0, 1e-4);
            TS_ASSERT_DELTA(node_forces[node_index][1], -5.0, 1e-4);
        }

        /* The last block of code provides an archiving test for the force class,
         * in a similar way to previous cell-based Chaste tutorials:
         *
         * Note that it is important to test archiving by using an abstract
         * pointer, so that you check that boost can identify and record which
         * concrete class it should be dealing with.
         * This tests the CHASTE_CLASS_EXPORT(MyForce) lines are implemented correctly.
         */
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "my_force.arch";
        {
            AbstractForce<2>* const p_force = new MyForce(2.6);
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_force;
            delete p_force;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractForce<2>* p_force;
            input_arch >> p_force;

            TS_ASSERT_DELTA(dynamic_cast<MyForce*>(p_force)->GetStrength(), 2.6, 1e-4);

            delete p_force;
        }

        /* We conclude the test by calling Destroy() on any singleton classes. */
        SimulationTime::Destroy();
    }

    /*
     * == Using the force in a cell-based simulation ==
     *
     * We now provide a test demonstrating how {{{MyForce}}} can be used
     * in a cell-based simulation.
     */
    void TestOffLatticeSimulationWithMyBoundaryCondition() throw(Exception)
    {
        /* The first thing to do, as before, is to set up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Once again we create a {{{MeshBasedCellPopulation}}}. */
        HoneycombMeshGenerator generator(20, 20, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We then pass in the cell population into a {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMyForce");
        simulator.SetEndTime(1.0);

        /* We create our force law and pass it to the {{{OffLatticeSimulation}}}. */
        MyForce force(0.5);
        simulator.AddForce(&force);

        /* We test that the Solve() method does not throw any exceptions. */
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        /* We conclude the test by calling {{{Destroy()}}} on any singleton classes. */
        SimulationTime::Destroy();
    }
};

#endif /*TESTCREATINGANDUSINGANEWFORCETUTORIAL_HPP_*/
