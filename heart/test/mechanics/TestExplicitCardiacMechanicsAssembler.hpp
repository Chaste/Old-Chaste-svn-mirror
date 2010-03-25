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


#ifndef TESTEXPLICITCARDIACMECHANICSASSEMBLER_HPP_
#define TESTEXPLICITCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "ExplicitCardiacMechanicsAssembler.hpp"
#include "ImplicitCardiacMechanicsAssembler.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraturePointsGroup.hpp"
#include "NonlinearElasticityTools.hpp"
#include "ReplicatableVector.hpp"

class TestExplicitCardiacMechanicsAssembler : public CxxTest::TestSuite
{
public:
    void TestWithSimpleContractionModel() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 4, 4);
        MooneyRivlinMaterialLaw<2> law(1);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        // NONPHYSIOL1 => NonphysiologicalContractionModel 1
        ExplicitCardiacMechanicsAssembler<2> assembler(NONPHYSIOL1,&mesh,"TestExplicitCardiacMech",fixed_nodes,&law);

        // coverage
        QuadraturePointsGroup<2> quad_points(mesh, *(assembler.GetQuadratureRule()));

        std::vector<double> calcium_conc(assembler.GetTotalNumQuadPoints(), 0.0);
        std::vector<double> voltages(assembler.GetTotalNumQuadPoints(), 0.0);

        assembler.SetCalciumAndVoltage(calcium_conc, voltages);

        // solve UP TO t=0. So Ta(lam_n,t_{n+1})=5*sin(0)=0, ie no deformation
        assembler.Solve(-0.01,0.0,0.01);
        TS_ASSERT_EQUALS(assembler.GetNumNewtonIterations(),0u);

        assembler.Solve(0.24,0.25,0.01);

        TS_ASSERT_DELTA(assembler.rGetDeformedPosition()[4](0),  0.8730, 1e-2);
        TS_ASSERT_DELTA(assembler.rGetDeformedPosition()[4](1), -0.0867, 1e-2);
    }

    // with stretch (and stretch-rate) independent contraction models the implicit and explicit schemes
    // are identical
    void TestCompareImplicitAndExplicitWithStretchIndependentContractionModel() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 4, 4);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        // NONPHYSIOL 1 - contraction model is of the form sin(t)
        ExplicitCardiacMechanicsAssembler<2> expl_solver(NONPHYSIOL1,&mesh,""/*"TestCompareExplAndImplCardiacSolvers_Exp"*/,fixed_nodes,&law);
        ImplicitCardiacMechanicsAssembler<2> impl_solver(NONPHYSIOL1,&mesh,""/*"TestCompareExplAndImplCardiacSolvers_Imp"*/,fixed_nodes,&law);

        double dt = 0.25;
        for(double t=0; t<3; t+=dt)
        {
            expl_solver.Solve(t,t+dt,dt);
            impl_solver.Solve(t,t+dt,dt);

            // computations should be identical
            TS_ASSERT_EQUALS(expl_solver.GetNumNewtonIterations(), impl_solver.GetNumNewtonIterations());
            for(unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](0),  impl_solver.rGetDeformedPosition()[i](0), 1e-10);
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](1),  impl_solver.rGetDeformedPosition()[i](1), 1e-10);
            }
        }
    }


    // with stretch-dependent contraction models the implicit and explicit schemes can be similar
    void TestCompareImplicitAndExplicitWithStretchDependentContractionModel() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 4, 4);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        // NONPHYSIOL 2 - contraction model is of the form lam*sin(t)
        ExplicitCardiacMechanicsAssembler<2> expl_solver(NONPHYSIOL2,&mesh,"TestCompareExplAndImplCardiacSolversStretch_Exp",fixed_nodes,&law);
        ImplicitCardiacMechanicsAssembler<2> impl_solver(NONPHYSIOL2,&mesh,"TestCompareExplAndImplCardiacSolversStretch_Imp",fixed_nodes,&law);

        expl_solver.WriteOutput(0);
        impl_solver.WriteOutput(0);

        unsigned counter = 1;

        double t0 = 0.0;
        double t1 = 0.25; // to be quite quick (min stretch ~=0.88), make this 5/4 (?) say for min stretch < 0.7
        double dt = 0.025;

        for(double t=t0; t<t1; t+=dt)
        {
            //std::cout << "\n **** t = " << t << " ****\n" << std::flush;

            expl_solver.SetWriteOutput(false);
            expl_solver.Solve(t,t+dt,dt);
            expl_solver.SetWriteOutput();
            expl_solver.WriteOutput(counter);

            impl_solver.SetWriteOutput(false);
            impl_solver.Solve(t,t+dt,dt);
            impl_solver.SetWriteOutput();
            impl_solver.WriteOutput(counter);

            // the solutions turn out to be very close to each other
            for(unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](0),  impl_solver.rGetDeformedPosition()[i](0), 2e-3);
                TS_ASSERT_DELTA(expl_solver.rGetDeformedPosition()[i](1),  impl_solver.rGetDeformedPosition()[i](1), 2e-3);
            }

            // visualisation in matlab or octave
            // run from CHASTETESTOUTPUT
            //
            // close all; figure; hold on
            // for i=0:10, x1 = load(['TestCompareExplAndImplCardiacSolversStretch_Exp/solution_',num2str(i),'.nodes']); plot(i,x1(5,1),'b*'); end
            // for i=0:10, x2 = load(['TestCompareExplAndImplCardiacSolversStretch_Imp/solution_',num2str(i),'.nodes']); plot(i,x2(5,1),'r*'); end
            //
            // close all; figure; hold on
            // for i=0:10, x2 = load(['TestCompareExplAndImplCardiacSolversStretch_Imp/solution_',num2str(i),'.nodes']); plot(x2(:,1),x2(:,2),'r*'); end

            counter++;
        }

        // check there was significant deformation - node 4 is (1,0)
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_LESS_THAN(impl_solver.rGetDeformedPosition()[4](0), 0.9);
    }

    // cover all other contraction model options which are allowed but not been used in a test so far
    void TestCoverage() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1, 1);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ExplicitCardiacMechanicsAssembler<2> expl_solver(NONPHYSIOL3,&mesh,"",fixed_nodes,&law);

        ExplicitCardiacMechanicsAssembler<2> expl_solver_with_nash(NASH2004,&mesh,"",fixed_nodes,&law);
        ExplicitCardiacMechanicsAssembler<2> expl_solver_with_kerchoffs(KERCHOFFS2003,&mesh,"",fixed_nodes,&law);

        // bad contraction model
        TS_ASSERT_THROWS_THIS(ExplicitCardiacMechanicsAssembler<2> assembler(NHS,&mesh,"",fixed_nodes,&law), "Unknown or stretch-rate-dependent contraction model");
    }
};

#endif /*TESTEXPLICITCARDIACMECHANICSASSEMBLER_HPP_*/
