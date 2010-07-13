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


#ifndef TESTCARDIACELECTROMECHANICSPROBLEM_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "NumericFileComparison.hpp"
#include "Hdf5DataReader.hpp"

class TestCardiacElectroMechanicsProblem : public CxxTest::TestSuite
{
public:
    void TestDeterminingWatchedNodes() throw(Exception)
    {
        HeartEventHandler::Disable();

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(NHS,
                                                     1.0, /* width (cm) */
                                                     1,   /* mech elem each dir */
                                                     96, /* elec elem each dir */
                                                     &cell_factory,
                                                     1.0, /* end time */
                                                     0.01, /* electrics timestep (ms) */
                                                     100, /* 100*0.01ms mech dt */
                                                     0.01,/* contraction model ode timestep */
                                                     "");

        c_vector<double,2> pos;
        pos(0) = 1.0;
        pos(1) = 0.0;
        problem.SetWatchedPosition(pos);
        problem.SetNoElectricsOutput();
        problem.Initialise();

        // have checked these hardcoded values correspond to the nodes
        // at (1,0);
        TS_ASSERT_EQUALS(problem.mWatchedElectricsNodeIndex, 96u);
        TS_ASSERT_EQUALS(problem.mWatchedMechanicsNodeIndex, 1u);

        //// would like to do the following....
        //CardiacElectroMechanicsProblem<2> problem2(&cell_factory,
        //                                           1, 10, 100, 0.01,
        //                                           "");
        //pos(1) = 1.1;
        //problem2.SetWatchedPosition(pos);
        //TS_ASSERT_THROWS_THIS(problem2.Initialise(), "");
        //// ... but the exception causes a segmentation fault and had to be replaced
        //// with an assert(0);
    }


    void TestImplicitNhs2dOneMechanicsElement() throw(Exception)
    {
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(NHS,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     10.0, /* end time */
                                                     0.01, /* electrics timestep (ms) */
                                                     100,  /* 100*0.01ms mech dt */
                                                     0.01, /* contraction model ode timestep */
                                                     "TestCardiacElectroMechOneElement");
        c_vector<double,2> pos;
        pos(0) = 0.05;
        pos(1) = 0.0;

        problem.SetWatchedPosition(pos);
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[1](0), 0.0497, 1e-4);

        OutputFileHandler handler("TestCardiacElectroMechOneElement",false);

        NumericFileComparison comparer(handler.GetOutputDirectoryFullPath() + "watched.txt","heart/test/data/good_watched.txt");
        TS_ASSERT(comparer.CompareFiles(1e-2));

        // check electrics output was written
        std::string command = "ls " + handler.GetOutputDirectoryFullPath() + "/electrics";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        // coverage
        CardiacElectroMechProbRegularGeom<2> prob_with_bad_model(NONPHYSIOL1,0.05,1,5,&cell_factory,1,0.01,100,0.01,"");
        TS_ASSERT_THROWS_CONTAINS(prob_with_bad_model.Solve(),"Invalid");

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    void TestWithKerchoffs() throw(Exception)
    {
        HeartEventHandler::Disable();

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(KERCHOFFS2003,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     20,    /* end time */     // dies at 7.28 with explicit (now using implicit)
                                                     0.01, /* electrics timestep (ms) */
                                                     100,   /* n times 0.01ms mech dt */
                                                     0.01, /* Kerchoffs ode timestep */
                                                     "TestCardiacEmWithKerchoffs");

        c_vector<double,2> pos;
        pos(0) = 0.05;
        pos(1) = 0.0;

        problem.SetWatchedPosition(pos);
        problem.SetNoElectricsOutput();
        problem.Initialise();

        problem.Solve();

        //visualise to verify

        // hardcoded result
        TS_ASSERT_EQUALS(problem.mWatchedMechanicsNodeIndex, 1u);
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](0), 0.0479, 0.0002);
    }

    //
    //  BAD test - fails with HYPRE (for some reason HYPRE can't solve the one of the linear systems, and 
    //  the search direction in the end doesn't decrease the residual), and also with ILU if you increase
    //  the number of elements (whether LR91 or N98 is used). Probably the active tension is too high. 
    //
    //
    void TestExplicitSolverWithNash2004() throw(Exception)
    {
        HeartEventHandler::Disable();

        PlaneStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(NASH2004,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     20,    /* end time */
                                                     0.01, /* electrics timestep (ms) */
                                                     100,   /* n times 0.01ms mech dt */
                                                     0.01, /*  ode timestep */
                                                     "TestExplicitWithNash");

        // coverage, this file is just X-direction fibres
        problem.SetVariableFibreSheetDirectionsFile("heart/test/data/1by1mesh_fibres.ortho", false);

        c_vector<double,2> pos;
        pos(0) = 0.05;
        pos(1) = 0.0;

        problem.SetWatchedPosition(pos);
        problem.SetNoElectricsOutput();
        problem.Initialise();

        problem.Solve();

        //visualise to verify

        // hardcoded result
        TS_ASSERT_EQUALS(problem.mWatchedMechanicsNodeIndex, 1u);
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](0), 0.0419, 0.0002);
    }
    
    // Sets up a short simulation on a square with zero stimulus, but a model with stretch activated channels.
    // Hacks the mechanics initial condition to correspond to some stretch, which should create a bit of 
    // SAC activity and increased voltage
    void TestWithMechanoElectricFeedback() throw (Exception)
    {
        PlaneStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 2> cell_factory(0.0);

        // set up two meshes of 1mm by 1mm by 1mm
        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.02, 0.1, 0.1);
 
        QuadraticMesh<2> mechanics_mesh(0.1, 0.1, 1, 1);

        // fix the nodes on x=0
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh,0,0);

        CardiacElectroMechanicsProblem<2> problem(NASH2004,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  fixed_nodes,
                                                  &cell_factory,
                                                  1,   /* end time */
                                                  0.01, /* electrics timestep (ms) */
                                                  100,  /* 100*0.01ms mech dt */
                                                  1.0,  /* contraction model ode dt */
                                                  "TestNobleSacActivatedByStretchTissue");

        // use MEF
        problem.UseMechanoElectricFeedback();

        problem.Initialise();

        // hack into the mechanics assembler and set up the current solution so that it corresponds to
        // the square of tissue being stretched
        //
        // Note after one timestep the tissue will have returned to the resting state as there are no
        // forces and no way at the moment of passing fixed-displacement boundary conditions down to the mech
        // assembler. 
        for(unsigned i=0; i<problem.mpMechanicsMesh->GetNumNodes(); i++)
        {
            double X = problem.mpMechanicsMesh->GetNode(i)->rGetLocation()[0];
            double Y = problem.mpMechanicsMesh->GetNode(i)->rGetLocation()[1];
            problem.mpCardiacMechAssembler->rGetCurrentSolution()[2*i]   = X*0.2;
            problem.mpCardiacMechAssembler->rGetCurrentSolution()[2*i+1] = Y*(1.0/1.2 - 1);
        }
            
                
        problem.Solve();
        
        // Get the voltage at the start and end of the simulation, check the stretch was passed down to the 
        // cell model and caused increased voltage
        
        Hdf5DataReader reader("TestNobleSacActivatedByStretchTissue/electrics", "voltage");
        Vec start_voltage = PetscTools::CreateVec(36);
        Vec end_voltage = PetscTools::CreateVec(36);
        reader.GetVariableOverNodes(start_voltage, "V", 0);
        reader.GetVariableOverNodes(end_voltage, "V", 1);
        ReplicatableVector start_voltage_repl(start_voltage);
        ReplicatableVector end_voltage_repl(end_voltage);

        for(unsigned i=0; i<start_voltage_repl.GetSize(); i++)
        {
            TS_ASSERT_LESS_THAN(start_voltage_repl[i], -90.0);
            TS_ASSERT_LESS_THAN(-90, end_voltage_repl[i]);
        }
        
        VecDestroy(start_voltage);         
        VecDestroy(end_voltage);         
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
