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


#ifndef TESTCARDIACELECTROMECHANICSPROBLEM_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "NumericFileComparison.hpp"
#include "Hdf5DataReader.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "MooneyRivlinMaterialLaw.hpp"

class TestCardiacElectroMechanicsProblem : public CxxTest::TestSuite
{
public:
    void TestDeterminingWatchedNodes() throw(Exception)
    {
        HeartEventHandler::Disable();

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     NHS,
                                                     1.0, /* width (cm) */
                                                     1,   /* mech elem each dir */
                                                     96,  /* elec elem each dir */
                                                     &cell_factory,
                                                     1.0,  /* end time */
                                                     0.01, /* electrics timestep (ms) */
                                                     1.0,  /* mechanics solve timestep */
                                                     0.01, /* contraction model ode timestep */
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
        
        // coverage
        TS_ASSERT_EQUALS( problem.GetSolidMechanicsProblemDefinition()->rGetFixedNodes()[0], 0u);
    }


    void TestImplicitNhs2dOneMechanicsElement() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     NHS,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     10.0, /* end time */
                                                     0.01, /* electrics timestep (ms) */
                                                     1.0,  /* mechanics solve timestep */
                                                     0.01, /* contraction model ode timestep */
                                                     "TestCardiacElectroMechOneElement");
        c_vector<double,2> pos;
        pos(0) = 0.05;
        pos(1) = 0.0;

        // cover SetMaterialLaw() - pass in the law that would be used anyway.
        NashHunterPoleZeroLaw<2> law;
        problem.SetMaterialLaw(&law);

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
        CardiacElectroMechProbRegularGeom<2> prob_with_bad_model(INCOMPRESSIBLE,NONPHYSIOL1,0.05,1,5,&cell_factory,1,0.01,1,0.01,"");
        TS_ASSERT_THROWS_CONTAINS(prob_with_bad_model.Solve(),"Invalid contraction model");

        TS_ASSERT_THROWS_CONTAINS(CardiacElectroMechProbRegularGeom<2> prob_with_bad_model(INCOMPRESSIBLE,NHS,0.05,1,5,&cell_factory,1,0.01,0.025,0.01,""),"does not divide");


        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    void TestWithKerchoffs() throw(Exception)
    {
        HeartEventHandler::Disable();

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     KERCHOFFS2003,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     20,    /* end time */     // dies at 7.28 with explicit (now using implicit)
                                                     0.01,  /* electrics timestep (ms) */
                                                     1.0,   /* mechanics solve timestep */
                                                     0.01,  /* Kerchoffs ode timestep */
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
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](1),-0.0003, 0.0002);
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

#ifdef MECH_USE_HYPRE
        TS_FAIL("This test is known to fail with HYPRE - see comments in test");
        return;
#endif

        PlaneStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     NASH2004,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     20,    /* end time */
                                                     0.01,  /* electrics timestep (ms) */
                                                     1.0,   /* mechanics solve timestep */
                                                     0.01,  /* nash ode timestep */
                                                     "TestExplicitWithNash");

        // coverage, this file is just X-direction fibres
        problem.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/1by1mesh_fibres.ortho", false);

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
 
        QuadraticMesh<2> mechanics_mesh(0.1, 0.1, 0.1);

        // fix the nodes on x=0
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh,0,0);

        CardiacElectroMechanicsProblem<2> problem(INCOMPRESSIBLE,
                                                  NASH2004,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  fixed_nodes,
                                                  &cell_factory,
                                                  1,    /* end time */
                                                  0.01, /* electrics timestep (ms) */
                                                  1.0,  /* mechanics solve timestep */
                                                  1.0,  /* contraction model ode dt */
                                                  "TestNobleSacActivatedByStretchTissue");

        // use MEF
        problem.UseMechanoElectricFeedback();

        problem.Initialise();

        // hack into the mechanics solver and set up the current solution so that it corresponds to
        // the square of tissue being stretched
        for(unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
        {
            double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
            double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
            problem.mpMechanicsSolver->rGetCurrentSolution()[2*i]   = X*0.2;
            problem.mpMechanicsSolver->rGetCurrentSolution()[2*i+1] = Y*(1.0/1.2 - 1);
        }

        // we are going to get the modified conductivity tensor directly, without (initially) calling solve,
        // so need to do the following, which is normally done inside the Solve
        problem.mpCardiacMechSolver->ComputeDeformationGradientAndStretchInEachElement(problem.mDeformationGradientsForEachMechanicsElement, problem.mStretchesForEachMechanicsElement);

        // just get the default conductivity so don't to hardcode it (1.75 at the moment)
        c_vector<double, 2> conductivities;
        HeartConfig::Instance()->GetIntracellularConductivities(conductivities);
        assert((conductivities(0)-conductivities(1))<1e-8);
        double default_conductivity = conductivities(0);


        // test directly that the conductivity is being computed using the deformation
        for(unsigned i=0; i<electrics_mesh.GetNumElements(); i++)
        {
            // sigma = F^{-1} sigma_undef F^{-T}, F=diag(1.2, 1.0/1.2), sigma = diag(1.75,1.75).
            const c_matrix<double,2,2>& r_tensor = problem.mpMonodomainProblem->GetMonodomainTissue()->rGetIntracellularConductivityTensor(i);
            TS_ASSERT_DELTA(r_tensor(0,0), default_conductivity/(1.2*1.2), 1e-9);
            TS_ASSERT_DELTA(r_tensor(0,1), 0.0,                            1e-9);
            TS_ASSERT_DELTA(r_tensor(1,0), 0.0,                            1e-9);
            TS_ASSERT_DELTA(r_tensor(1,1), default_conductivity*(1.2*1.2), 1e-9);
        }

                
        // Note after one timestep the tissue will have returned to the resting state as there are no
        // forces and no way at the moment of passing fixed-displacement boundary conditions down to the mech
        // solver.
        problem.Solve();
        
        // Get the voltage at the start and end of the simulation, check the stretch was passed down to the 
        // cell model and caused increased voltage
        
        Hdf5DataReader reader("TestNobleSacActivatedByStretchTissue/electrics", "voltage");
        unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();
        TS_ASSERT_EQUALS(num_timesteps, 2u);
        
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

    // Similar to first part of above test, except the deformation isn't just constant stretch here, it is
    // different in the two elements of the mechanics mesh
    void TestWithMechanoElectricFeedbackHeterogeneousStretch() throw (Exception)
    {
        // irrelevant, not going to call solve
        PlaneStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 2> cell_factory(0.0);

        // set up two meshes of 1mm by 1mm by 1mm
        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.025, 0.1, 0.1); // the h should be such that there are an even number of elements per dim, so
                                                                  // electrical element centroids don't lie on mechanics element boundaries
                                                                  // (for the below test to not have to worry about boundary cases).

        QuadraticMesh<2> mechanics_mesh(0.1, 0.1, 0.1); // 2 elements

        // irrelevant, not going to call solve
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh,0,0);

        CardiacElectroMechanicsProblem<2> problem(INCOMPRESSIBLE,
                                                  NASH2004,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  fixed_nodes,
                                                  &cell_factory,
                                                  1,    /* end time */
                                                  0.01, /* electrics timestep (ms) */
                                                  1.0,  /* mechanics solve timestep */
                                                  1.0,  /* contraction model ode dt */
                                                  "");

        // use MEF
        problem.UseMechanoElectricFeedback();

        problem.Initialise();

        // hack into the mechanics solver and set up the current solution so that it corresponds to
        // the some stretch in the upper element
        for(unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
        {
            double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
            double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
            problem.mpMechanicsSolver->rGetCurrentSolution()[2*i]   = X*Y*0.1; // displacement is zero except for (1,1) node
            problem.mpMechanicsSolver->rGetCurrentSolution()[2*i+1] = X*Y*0.1; // displacement is zero except for (1,1) node
        }

        // we are going to get the modified conductivity tensor directly, without (initially) calling solve,
        // so need to do the following, which is normally done inside the Solve
        problem.mpCardiacMechSolver->ComputeDeformationGradientAndStretchInEachElement(problem.mDeformationGradientsForEachMechanicsElement, problem.mStretchesForEachMechanicsElement);

        c_matrix<double,2,2> F;
        F(0,0) = 1.01;
        F(1,0) = 0.01;
        F(0,1) = 0.01;
        F(1,1) = 1.01;
        c_matrix<double,2,2> inverse_C = prod(Inverse(F), trans(Inverse(F)));

        // just get the default conductivity so don't to hardcode it (1.75 at the moment)
        c_vector<double, 2> conductivities;
        HeartConfig::Instance()->GetIntracellularConductivities(conductivities);
        assert((conductivities(0)-conductivities(1))<1e-8);
        double default_conductivity = conductivities(0);

        // test directly that the conductivity is being computed using the deformation
        for(unsigned i=0; i<electrics_mesh.GetNumElements(); i++)
        {
            // sigma = F^{-1} sigma_undef F^{-T},
            const c_matrix<double,2,2>& r_tensor = problem.mpMonodomainProblem->GetMonodomainTissue()->rGetIntracellularConductivityTensor(i);
            c_vector<double,2> centroid = electrics_mesh.GetElement(i)->CalculateCentroid();
            if(centroid(0)+centroid(1)<0.1)
            {
                // in mechanics element with corners (0,0),(0,1),(1,0) -- F=I here
                TS_ASSERT_DELTA(r_tensor(0,0), default_conductivity, 1e-9);
                TS_ASSERT_DELTA(r_tensor(0,1), 0.0,                  1e-9);
                TS_ASSERT_DELTA(r_tensor(1,0), 0.0,                  1e-9);
                TS_ASSERT_DELTA(r_tensor(1,1), default_conductivity, 1e-9);
            }
            else
            {
                // in mechanics element with corners (0,1),(1,0),(1,1) -- F \ne I here.
                // sigma_def = F^{-1} sigma_undef F^{-T}, but this is 1.75 C^{-1} since sigma_undef = 1.75*I
                TS_ASSERT_DELTA(r_tensor(0,0), default_conductivity*inverse_C(0,0), 1e-9);
                TS_ASSERT_DELTA(r_tensor(0,1), default_conductivity*inverse_C(0,1), 1e-9);
                TS_ASSERT_DELTA(r_tensor(1,0), default_conductivity*inverse_C(1,0), 1e-9);
                TS_ASSERT_DELTA(r_tensor(1,1), default_conductivity*inverse_C(1,1), 1e-9);
            }
        }
    }

    void TestWithCompressibleApproach() throw(Exception)
    {
        EXIT_IF_PARALLEL; // #1913 currently, the compressible preconditioner is ICC, which is only supported in sequential

        HeartEventHandler::Disable();

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(COMPRESSIBLE,
                                                     KERCHOFFS2003,
                                                     0.05, /* width (cm) */
                                                     1,    /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     20,    /* end time */
                                                     0.01,  /* electrics timestep (ms) */
                                                     1.0,   /* mechanics solve timestep */
                                                     0.01,  /* Kerchoffs ode timestep */
                                                     "TestCompressibleWithKerchoffs");

        problem.Solve();

        // Mainly just testing no errors when Solve was called.
        // The results of this test can be visually compared with the results of the
        // equivalent incompressible simulation in TestWithKerchoffs.

        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](0), 0.0465, 0.0002);
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[1](1),-0.0012, 0.0002);
    }

    void TestCardiacElectroMechanicsHeterogeneousMaterialLaws() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-5000*1000);

        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.02/*stepsize*/, 0.1/*length*/, 0.1/*width*/);

        QuadraticMesh<2> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.02, 0.1, 0.1 /*as above with a different stepsize*/);

        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 0, 0.0);

        CardiacElectroMechanicsProblem<2> problem(INCOMPRESSIBLE,
                                                  KERCHOFFS2003,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  fixed_nodes,
                                                  &cell_factory,
                                                  20,   // end time
                                                  0.01, // electrics timestep (ms)
                                                  1.0,  // mechanics solve timestep
                                                  0.01,  // contraction model ode timestep
                                                  "TestCardiacElectroMechanicsHeterogeneousMaterialLaws" /* output directory */);


        /* Create two laws to have stiff and soft tissue */
        std::vector<AbstractMaterialLaw<2>*> law;
        MooneyRivlinMaterialLaw<2> stiff_law(1.0);
        MooneyRivlinMaterialLaw<2> soft_law(1.0/5.0);
        for (TetrahedralMesh<2,2>::ElementIterator iter = mechanics_mesh.GetElementIteratorBegin();
             iter != mechanics_mesh.GetElementIteratorEnd();
             ++iter)
        {
            if (((iter)->CalculateCentroid()[1] >= 0.04)
                 && ((iter)->CalculateCentroid()[1] <= 0.06))
            {
                law.push_back(&stiff_law);
            }
            else
            {
                law.push_back(&soft_law);
            }
        }

        problem.SetMaterialLaw(law);
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();

        // node 5 starts at (1,0)
        assert(fabs(mechanics_mesh.GetNode(5)->rGetLocation()[0] - 0.1)<1e-6);
        assert(fabs(mechanics_mesh.GetNode(5)->rGetLocation()[1])<1e-6);
///\todo #1948
        // Visualised solution to check heterogeneous stiffnesses are taken into account,
        // here we just have a hardcoded test to check nothing has changed
        TS_ASSERT_DELTA(r_deformed_position[5](0),  0.0910, 1e-4);
        TS_ASSERT_DELTA(r_deformed_position[5](1), -0.0039, 1e-4);
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
