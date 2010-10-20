

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
#ifndef TESTCARDIACCHECKPOINTINGANDRESTARTINGTUTORIAL_HPP_
#define TESTCARDIACCHECKPOINTINGANDRESTARTINGTUTORIAL_HPP_
/*
 * = Checkpointing and restarting cardiac simulations =
 * 
 * In this tutorial we show how to save and reload cardiac simulations
 *
 * EMPTYLINE
 * 
 * `CardiacSimulationArchiver` has to be included. Archiving includes often have to be included first.
 */
#include <cxxtest/TestSuite.h>
#include "CardiacSimulationArchiver.hpp"
#include "BidomainProblem.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PlaneStimulusCellFactory.hpp"


class TestCardiacCheckpointingAndRestartingTutorial : public CxxTest::TestSuite
{
public:
    /* First, the checkpointing test. */
    void TestCheckpointing() throw(Exception)
    {
        /* We set up exactly the same simulation as in UserTutorials/AnotherBidomainSimulation */
        HeartConfig::Instance()->Reset();
        
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML,2> cell_factory(-2000000);
        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainCheckpointingTutorial");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_800_elements", cp::media_type::Orthotropic);

        double scale = 2;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75*scale, 0.19*scale));
        BidomainProblem<2> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        /* To save the entire simulation, use the `CardiacSimulationArchiver` class, as shown in the following. 
         * Note the `BidomainProblem<2>` as the template parameter. The output directory is relative to 
         * CHASTE_TEST_OUTPUT. */
        CardiacSimulationArchiver<BidomainProblem<2> >::Save(bidomain_problem, "BidomainCheckpointingTutorial/saved_simulation"); 
    }


 
    /* This is how to restart the test. */
    void TestRestarting() throw(Exception)
    {
        /* To restart from the saved simulation directory we  use the `CardiacSimulationArchiver` class, as shown in the following. 
         * Note the `BidomainProblem<2>` as the template parameter again.  The dimension (2) must match the one given in the 
         * saved archive directory.  
         * The output directory is again relative to CHASTE_TEST_OUTPUT. */
        BidomainProblem<2>* p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<2> >::Load("BidomainCheckpointingTutorial/saved_simulation");

        /* The simulation duration has to be amended.
         * Note that the duration is always given with respect to the origin of the first solve.
         * This means that we are running from {{{t=5 ms}}} (the end of the previous simulation) to {{{t=10 ms}}}.
         * The output files are concatenated so that they appear to be made by a single simulation running from
         * {{{t=0 ms}}} to {{{t=10 ms}}}.
         */
        HeartConfig::Instance()->SetSimulationDuration(10); //ms
        
        /* The point of checkpointing and restarting is that there may be something which we want to change
         * during the course of experiment.  Here we change the conductivity. */
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(3.0, 0.3));

        p_bidomain_problem->Solve();
    }
};
     



#endif /*TESTCARDIACCHECKPOINTINGANDRESTARTINGTUTORIAL_HPP_*/
