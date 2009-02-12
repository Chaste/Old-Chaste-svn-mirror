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
#ifndef TESTRUNNINGBIDOMAINSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGBIDOMAINSIMULATIONSTUTORIAL_HPP_
/*
 * = An example showing how to run bidomain simulations (for monodomain, it is essentially the same) =
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
 * implement the method {{{CreateCardiacCellForTissueNode(unsigned nodeNum)}}}, which returns
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
 * class will set this before it calls {{{CreateCardiacCellForTissueNode}}}, so it can be used
 * in that method); {{{mTimestep}}}, a double (see below); and {{{mpSolver}}} a forward
 * euler ode solver (see below). */
private:
    SimpleStimulus *mpStimulus;

public:
    /* Our contructor takes in nothing. It calls the constructor of
     * {{{AbstractCardiacCellFactory}}} with 0.01 - this is what {{{mTimestep}}} will be set
     * to. We also initialise the stimulus to have magnitude -6000 '''TODO:units''' and duration 0.5ms.
     */
    PointStimulus2dCellFactory() : AbstractCardiacCellFactory<2>()
    {
        mpStimulus = new SimpleStimulus(-6000.0, 0.5);
    }

    /* Now we implement the pure method which needs to be implemented. We return
     * a LR91 cell for each node, with the node at (0,0) given the non-zero stimulus,
     * and all other nodes given the zero stimulus. Note that we use {{{mpMesh}}},
     * {{{mTimestep}}}, {{{mpZeroStimulus}}} and {{{mpSolver}}} which are all
     * members of the base class. The timestep and solver are defined in the base
     * class just so that the user doesn't have to create them here. */
    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->mpMesh->GetNode(nodeIndex)->rGetLocation()[0];
        double y = this->mpMesh->GetNode(nodeIndex)->rGetLocation()[1];
        if (fabs(x)+fabs(y)<1e-6) // ie if (x,y)==(0,0). An alternative would be if(norm_2(this->mpMesh->GetNode(nodeIndex)->rGetLocation())<1e-6)
        {
            /* Even if running a bidomain simulation, only the intra-cellular stimulus
             * should be given here.  There is a separate Electrodes class for applying
             * extra-cellular stimuli.
             */
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus);
        }
        else
        {
            /* The other cells have zero stimuli. */
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
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
        /* The {{{HeartConfig}}} class is used to set various parameters. It gets the default values 
         * from !ChasteDefaults.xml (except the values in the 'Simulation' block of the XML file,
         * which is only used by the Chaste executable). Parameters in this file can be re-set 
         * with {{{HeartConfig}}} if the user wishes, and other paramters such as end time must be set
         * using {{{HeartConfig}}}. Let us begin by setting the end time (in ms), the mesh to use, and the
         * output directory and filename-prefix.
         */ 
        HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/square_128_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainTutorial");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        
        /* Next, we have to create a cell factory of the type we defined above. */
        PointStimulus2dCellFactory cell_factory;

        /* Now we create a problem class using (a pointer to) the cell factory. */
        BidomainProblem<2> bidomain_problem( &cell_factory );
 
        /* This is enough setup to run a simulation: we could now call {{{Initialise()}}}
         * and {{{Solve()}}} to run... */
        // bidomain_problem.Initialise();
        // bidomain_problem.Solve();

        /* ..However, instead we show how to set a few more parameters. To set the conductivity values
         *  in the principal fibre, sheet and normal directions do the following.
         * Note that {{{Create_c_vector}}} is just a helper method for creating a {{{c_vector<double,DIM>}}}
         * of the correct size (2, in this case). Make sure these methods are called before 
         * {{{Initialise()}}}.
         */
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.19));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(6.2, 2.4));
        /* Let us also reset the surface-area-to-volume ratio and the capacitance */
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0); // TODO:units
        HeartConfig::Instance()->SetCapacitance(1.0); // TODO:units

        /* Now we call {{{Initialise()}}}... */
        bidomain_problem.Initialise();

        /* The output will be written to /tmp/USER_NAME/testoutput/BidomainTutorial 
         * in hdf5 format. If you want visualise the results afterwards, call the
         * following method now, and the mesh and output will be converted to meshalyzer 
         * format at the end of the simulation.
         */
        bidomain_problem.ConvertOutputToMeshalyzerFormat(); 

        /* Now we call Solve() to run the simulation.  
         * Note that if you want to view the progress of longer simulations
         * go to the the output directory and look at the file 
         * {{{progress_status.txt}}}, which will say the percentage of the 
         * simulation run. A useful linux command is therefore {{{watch tail progress_status.txt}}}
         * which will repeatedly display the last few lines of this file. */
        bidomain_problem.Solve();

        /* To now visualise the results, go to /tmp/USER_NAME/testoutput/BidomainTutorial/output,
         * where you should find the mesh and output, and run meshalyzer.
         * 
         * EMPTYLINE
         * 
         * The easiest way to look at the resultant voltage values (for the last timestep - 
         * the data for the previous timesteps is written to file but not retained) is to
         * use a {{{ReplicatableVector}}}. {{{bidomain_problem.GetSolution())}}} returns a !PetSc vector
         * of the form (V_0, phi_0, V_1, phi_e_1, ... V_n, phi_e_n), and we can create a 
         * {{{ReplicatableVector}}} for easy access to this !PetSc vector's data. (This won't be very
         * efficient with huge problems in parallel).
         */
        ReplicatableVector res_repl(bidomain_problem.GetSolution());
        for(unsigned i=0; i<res_repl.size(); i++)
        {
        //    std::cout << res_repl[i] << "\n";
        }

        /* Alternatively, we show how to access the voltage values using the {{{DistributedVector}}}
         * class, which can be used to only iterate over the values of the voltage owned
         * by that process.
         */
        DistributedVector dist_bidomain_voltage(bidomain_problem.GetSolution());
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
    }
};

#endif /*TESTRUNNINGBIDOMAINSIMULATIONSTUTORIAL_HPP_*/
