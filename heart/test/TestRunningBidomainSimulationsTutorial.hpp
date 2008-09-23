/*

Copyright (C) University of Oxford, 2008

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
#ifndef TESTRUNNINGBIDOMAINSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGBIDOMAINSIMULATIONSTUTORIAL_HPP_
/*
 * == Introduction ==
 *
 * In this tutorial we show how Chaste is used to run a standard bidomain simulation.
 * Note that monodomain simulations are run very similarly.
 *
 * The first thing that needs to be done, when writing any Chaste test,
 * is to include the following header.
 */

#include <cxxtest/TestSuite.h>
/* The main class to be used for running bidomain simulations is {{{BidomainProblem}}}. */
#include "BidomainProblem.hpp"
/* The {{{PlaneStimulusCellFactory}}} is a useful class to include (see later). */
#include "PlaneStimulusCellFactory.hpp"
/* {{{LuoRudyIModel1991OdeSystem}}} is the cell model which will be used in this simulation.*/
#include "LuoRudyIModel1991OdeSystem.hpp"
/* All tests which run cardiac simulations (which use Petsc) should include
 * {{{PetscSetupAndFinalize.hpp}}}.  This class ensures that {{{PetscInitialise()}}}
 * is called with the appropriate arguments before any tests in the suite are run. */
#include "PetscSetupAndFinalize.hpp"

/* EMPTYLINE
 *
 * == Defining a cell factory ==
 *
 * EMPTYLINE
 *
 * All mono/bidomain simulations need a ''cell factory'' as input. This is a class
 * which tells the problem class what type of cardiac cells to create. The cell-factory
 * class has to inherit from {{{AbstractCardiacCellFactory<DIM>}}}, which means it must
 * implement the method {{{CreateCardiacCellForNode(unsigned nodeNum)}}}, which returns
 * a pointer to an {{{AbstractCardiacCell}}}. Note, some concrete cell factories have
 * been defined, such as the {{{PlaneStimulusCellFactory}}}, which could be used in the
 * simulation, but for completeness we create our own cell factory in this test. For
 * complicated problems with, say, heterogeneous cell types or particular stimuli, a
 * new cell factory will have to be defined by the user for their particular problem.
 *
 * EMPTYLINE
 *
 * This cell factory is a simple cell factory where every cell is a Luo-Rudy 91 cell,
 * and only the cell at position (0,0) is given a non-zero stimulus.
 */
class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
/* Declare pointer to an {{{SimpleStimulus}}} for the cell which is stimulated.
 * Note that {{{AbstractCardiacCellFactory}}} also has as protected members: {{{mpZeroStimulus}}}
 * of type {{{ZeroStimulus}}}; {{{mpMesh}}}, a pointer to the mesh used (the problem
 * class will set this before it calls {{{CreateCardiacCellForNode}}}, so it can be used
 * in that method); {{{mTimestep}}}, a double (see below); and {{{mpSolver}}} a forward
 * euler ode solver (see below). */
private:
    SimpleStimulus *mpStimulus;

public:
    /* Our contructor takes in nothing. It calls the constructor of
     * {{{AbstractCardiacCellFactory}}} with 0.01 - this is what {{{mTimestep}}} will be set
     * to. We also initialise the stimulus to have magnitude -6000 and duration 0.5ms.
     */
    PointStimulus2dCellFactory() : AbstractCardiacCellFactory<2>()
    {
        mpStimulus = new SimpleStimulus(-6000.0, 0.5);
    }

    /* Now we implement the pure method which needs to be implemented. We return
     * a LR91 cell for each node, with the node at (0,0) given the non-zero stimulus,
     * and all other nodes given the zero stimulus. Note that we use {{{mpMesh}}},
     * {{{mTimestep}}}, {{{mpZeroStimulus}}} and {{{mpSolver}}} which are all
     * members of the base class. The timestep and solver being defined in the base
     * class are just so that the user doesn't have to create them here. */
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned nodeIndex)
    {
        double x = this->mpMesh->GetNode(nodeIndex)->rGetLocation()[0];
        double y = this->mpMesh->GetNode(nodeIndex)->rGetLocation()[1];
        if (fabs(x)+fabs(y)<1e-6) // ie if (x,y)==(0,0)
        {
            /* '''Note:''' As this is going to be used in a bidomain simulation, two
             * stimuli, the intra- and extra-cellular stimuli need to be given. We give the
             * cell a non-zero intra and zero extra-cellular stimulus. */
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus, mpZeroStimulus);
            /* For monodomain simulations we could have just done */
            //return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            /* The other cells have zero stimuli for both the intra and extra-cellular spaces. */
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus, mpZeroStimulus);
        }
    }

    /* The destructor just deletes the memory for the stimulus. Note the the problem
     * class deals with deleting the cells. */
    ~PointStimulus2dCellFactory()
    {
        delete mpStimulus;
    }
};

/*
 * EMPTYLINE
 *
 * == Running the simulation ==
 *
 * EMPTYLINE
 *
 * Now we can define the test class, which must inherit from {{{CxxTest::TestSuite}}}
 * as usual. */
class TestRunningBidomainSimulationsTutorial : public CxxTest::TestSuite
{
/* Tests should be public... */
public:
    /* Define the test. Note the {{{throw(Exception)}}} - without this exception messages
     * might not get printed out. */
    void TestSimpleSimulation() throw(Exception)
    {
        /* '''TODO - fill in comments on using HeartConfig''' */
        HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
        
        /* Next, we have to create a cell factory of the type we defined above. */
        PointStimulus2dCellFactory cell_factory;

        /* Now we create a problem class using (a pointer to) the cell factory. */
        BidomainProblem<2> bidomain_problem( &cell_factory );

        /* Next, some things which have to be set: the mesh filename, and the end time
         * (in ms). */
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/square_128_elements");
 
        /* If we want output to be written we need to set the output directory and output
         * file prefix.
         */
        HeartConfig::Instance()->SetOutputDirectory("BidomainTutorial");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        /* If this was enough setup, we could then call {{{Initialise()}}}
         * and {{{Solve()}}} to run the simulation... */
        // bidomain_problem.Initialise();
        // bidomain_problem.Solve();
        /* ..Instead we show how to set a few parameters. To
         * set the conductivity ''values'' in the principal fibre, sheet and normal directions do the following.
         * Note that {{{Create_c_vector}}} is just a helper method for creating a {{{c_vector<double,DIM>}}}
         * of the correct size (2, in this case). Note that these methods need to be called before
         * {{{Initialise()}}} '''is this true?''' '''todo - fix this'''*/
        //bidomain_problem.SetIntracellularConductivities(Create_c_vector(1.75, 0.19));
        //bidomain_problem.SetExtracellularConductivities(Create_c_vector(6.2, 2.4));
        
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.19));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(6.2, 2.4));
        
        /* Now we call {{{Initialise()}}}... */
        bidomain_problem.Initialise();
        /* .. and set the surface-area-to-volume ratio and capicitance. 
         */
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        /* Now we call Solve() to run the simulation. Output will be written
         * to /tmp/USER_NAME/testoutput/BidomainTutorial in hdf5 format. '''todo: To visualise...'''*/
        bidomain_problem.Solve();

        /* Finally, we show how to access the voltage values (at the final timestep, the
         * data for previous timesteps is not retained), using the {{{DistributedVector}}}
         * class. The call {{{bidomain_problem.GetVoltage())}}} returns a Petsc vector
         * of the form (V_0, phi_0, V_1, phi_e_1, ... V_n, phi_e_n), and the {{{DistributedVector}}}
         * class can be used to get the values. */
        DistributedVector dist_bidomain_voltage(bidomain_problem.GetVoltage());
        DistributedVector::Stripe bidomain_voltage(dist_bidomain_voltage, 0);
        DistributedVector::Stripe extracellular_potential(dist_bidomain_voltage, 1);

        /* A loop over all the components owned by this process.. */
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            /* .. and a simple test, that the 'last' node was stimulated: */
            if (index.Global==bidomain_problem.rGetMesh().GetNumNodes()-1) // ie if the last node
            {
                TS_ASSERT_LESS_THAN(0, bidomain_voltage[index]);
            }
        }

        /* Recall that the {{{ReplicatableVector}}} class can also be used for easier access. */
        //ReplicatableVector res_repl(bidomain_problem.GetVoltage());
        //for(unsigned i=0; i<res_repl.size(); i++)
        //{
        //    std::cout << res_repl[i] << "\n";
        //}
    }
};

#endif /*TESTRUNNINGBIDOMAINSIMULATIONSTUTORIAL_HPP_*/
