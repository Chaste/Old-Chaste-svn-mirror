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


#ifndef TESTIMPLICITCARDIACMECHANICSASSEMBLER_HPP_
#define TESTIMPLICITCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "ImplicitCardiacMechanicsAssembler.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraturePointsGroup.hpp"
#include "NonlinearElasticityTools.hpp"
#include "ReplicatableVector.hpp"

class TestImplicitCardiacMechanicsAssembler : public CxxTest::TestSuite
{
public:
    void TestCompareJacobians() throw(Exception)
    {
        EXIT_IF_PARALLEL; // this test usually passes in ||, but sometimes fails

        QuadraticMesh<2> mesh(1.0, 1.0, 1, 1);
        MooneyRivlinMaterialLaw<2> law(0.02);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ImplicitCardiacMechanicsAssembler<2> assembler(NHS, &mesh,"",fixed_nodes,&law);

        std::vector<double> calcium_conc(assembler.GetTotalNumQuadPoints(), 0.0);
        std::vector<double> voltages(assembler.GetTotalNumQuadPoints(), 0.0);

        for(unsigned i=0; i<calcium_conc.size(); i++)
        {
            calcium_conc[i] = 0.05;
        }

        assembler.SetCalciumAndVoltage(calcium_conc, voltages);

        // NOTE: calling CompareJacobians below bypasses calling Solve(t0,t1,dt), hence the
        // time info will not be set. We therefore must explicitly set them here.
        assembler.mCurrentTime = 0.0;
        assembler.mNextTime = 0.01;
        assembler.mOdeTimestep = 0.01;


        ///////////////////////////////////////////////////////////////////
        // compute numerical jacobian and compare with analytic jacobian
        // (about u=0, p=p0)
        ///////////////////////////////////////////////////////////////////
        assembler.AssembleSystem(true, true);
        ReplicatableVector rhs_vec(assembler.mpLinearSystem->rGetRhsVector());
        unsigned num_dofs = rhs_vec.GetSize();
        double h = 1e-6;
        int lo, hi;
        MatGetOwnershipRange(assembler.mpLinearSystem->rGetLhsMatrix(), &lo, &hi);

        for(unsigned j=0; j<num_dofs; j++)
        {
            assembler.mCurrentSolution.clear();
            assembler.FormInitialGuess();
            assembler.mCurrentSolution[j] += h;

            assembler.AssembleSystem(true, false);

            ReplicatableVector perturbed_rhs( assembler.mpLinearSystem->rGetRhsVector() );

            for(unsigned i=0; i<num_dofs; i++)
            {
                if((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = assembler.mpLinearSystem->GetMatrixElement(i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec[i])/h;
                    if((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-2);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-4);
                    }
                }
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);

        // coverage - test default material law works ok
        ImplicitCardiacMechanicsAssembler<2> another_assembler(NHS, &mesh,"",fixed_nodes);
        c_matrix<double,2,2> F = zero_matrix<double>(2,2);
        F(0,0)=F(1,1)=1.1;
        double pressure = 1;
        c_matrix<double,2,2> S;
        another_assembler.mMaterialLaws[0]->Compute1stPiolaKirchoffStress(F,pressure,S);
        TS_ASSERT_DELTA(S(0,0), 1.5805, 1e-3);
    }

    // A test where we specify the 'resting' intracellular calcium concentration
    // for which the active tension should be zero, so should solve in 0 newton
    // iterations
    void TestWithZeroActiveTension() throw(Exception)
    {
        EXIT_IF_PARALLEL; // this test usually passes in ||, but sometimes fails

        QuadraticMesh<2> mesh(1.0, 1.0, 8, 8);
        MooneyRivlinMaterialLaw<2> law(0.02);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ImplicitCardiacMechanicsAssembler<2> assembler(NHS, &mesh,"ImplicitCardiacMech/ZeroActiveTension",fixed_nodes,&law);

        TS_ASSERT_EQUALS(assembler.GetTotalNumQuadPoints(), mesh.GetNumElements()*9u);

        // 0.0002 is the initial Ca conc in Lr91
        std::vector<double> calcium_conc(assembler.GetTotalNumQuadPoints(), 0.0002);
        std::vector<double> voltages(assembler.GetTotalNumQuadPoints(), 0.0);
        assembler.SetCalciumAndVoltage(calcium_conc, voltages);

        assembler.Solve(0,0.1,0.01);

        TS_ASSERT_EQUALS(assembler.GetNumNewtonIterations(), 0u);
    }


    // This test compares against the deprecated explicit CardiacMechanicsAssembler
    // in the dealii folder. The implicit and explicit assemblers could be compared
    // by solving with the implicit for given Ca, finding the Ta (at each node)
    // which corresponds to this Ca, then checking the implicit solution gives zero
    // residual in the explicit assembler with those Ta's.
    //
    // There is no Chaste explicit cardiac mechanics assembler, so can't do the same
    // here. Instead we do exactly the same implicit solve as in TestCompareWithExplicit
    // in dealii/test/TestImplicitCardiacMechanicsAssembler.hpp (see r4395 say),
    // and check the two implicit solvers agree
    void TestCompareWithDeadExplicitSolver() throw(Exception)
    {
        EXIT_IF_PARALLEL; // unlike above tests, this one doesn't pass in parallel (only written for sequential)

        // note 8 elements is assumed in the fixed nodes
        QuadraticMesh<2> mesh(1.0, 1.0, 8, 8);
        MooneyRivlinMaterialLaw<2> law(0.02);

        // fix all nodes on elements containing the origin (as was done in
        // dealii test)
        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        fixed_nodes.push_back(82);
        fixed_nodes.push_back(9);
        fixed_nodes.push_back(1);
        fixed_nodes.push_back(81);

        ImplicitCardiacMechanicsAssembler<2> assembler(NHS, &mesh,"ImplicitCardiacMech/CompareWithExplicit",fixed_nodes,&law);

        std::vector<double> calcium_conc(assembler.GetTotalNumQuadPoints(), 1); // unrealistically large Ca (but note random material law used)
        std::vector<double> voltages(assembler.GetTotalNumQuadPoints(), 0.0);
        assembler.SetCalciumAndVoltage(calcium_conc, voltages);

        assembler.Solve(0,0.01,0.01);

        // Visually compared results, they are identical to the dealii results
        // Hardcoded value for (1,1) node
        TS_ASSERT_DELTA(assembler.rGetDeformedPosition()[80](0),  0.98822 /*dealii*/, 5e-4);
        TS_ASSERT_DELTA(assembler.rGetDeformedPosition()[80](1),  1.01177 /*dealii*/, 3e-4);
        // Hardcoded value for (0,1) node
        TS_ASSERT_DELTA(assembler.rGetDeformedPosition()[72](0), -0.00465 /*dealii*/, 4e-4);
        TS_ASSERT_DELTA(assembler.rGetDeformedPosition()[72](1),  1.00666 /*dealii*/, 1e-4);
    }


    // Specifies a non-constant active tension and checks the lambda behaves
    // as it should do. Also has hardcoded tests
    void TestSpecifiedActiveTensionCompression() throw(Exception)
    {
        EXIT_IF_PARALLEL; // unlike above tests, this one doesn't pass in parallel (only written for sequential)

        // NOTE: test hardcoded for num_elem = 4
        QuadraticMesh<2> mesh(1.0, 1.0, 4, 4);
        MooneyRivlinMaterialLaw<2> law(0.02);

        // need to leave the mesh as unfixed as possible
        std::vector<unsigned> fixed_nodes(2);
        fixed_nodes[0] = 0;
        fixed_nodes[1] = 5;

        ImplicitCardiacMechanicsAssembler<2> assembler(NHS, &mesh,"ImplicitCardiacMech/SpecifiedCaCompression",fixed_nodes,&law);
        QuadraturePointsGroup<2> quad_points(mesh, *(assembler.GetQuadratureRule()));

        std::vector<double> calcium_conc(assembler.GetTotalNumQuadPoints());
        for(unsigned i=0; i<calcium_conc.size(); i++)
        {
            double Y = quad_points.Get(i)(1);
            // 0.0002 is the initial Ca conc in Lr91, 0.001 is the greatest Ca conc
            // value in one of the Lr91 TestIonicModel tests
            calcium_conc[i] = 0.0002 + 0.001*Y;
        }

        std::vector<double> voltages(assembler.GetTotalNumQuadPoints(), 0.0);
        assembler.SetCalciumAndVoltage(calcium_conc, voltages);

        // solve for quite a long time to get some deformation
        assembler.Solve(0,10,1);

        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 24 is the top-right corner node,
        TS_ASSERT_DELTA( assembler.rGetDeformedPosition()[24](0), 0.9413, 1e-3);
        TS_ASSERT_DELTA( assembler.rGetDeformedPosition()[24](1), 1.0582, 1e-3);

        std::vector<double>& lambda = assembler.rGetFibreStretches();

        // the lambdas should be less than 1 (ie compression), and also
        // should be near the same for any particular value of Y, ie the
        // same along any fibre. Lambda should decrease nonlinearly.
        // Uncomment trace and view in matlab (plot y against lambda)
        // to observe this - SEE FIGURE ATTACHED TO TICKET #757.
        // The lambda are constant for given Y if Y>0.1.5 (ie not near fixed nodes)
        // and a cubic polynomial can be fitted with matlab
        for(unsigned i=0; i<lambda.size(); i++)
        {
            TS_ASSERT_LESS_THAN(lambda[i], 1.0);

            // Get the value of Y for the point
            double Y = quad_points.Get(i)(1);
            // Lambda should be near a value obtained by fitting a
            // cubic polynomial of lambda against Y
            // Matlab code:
            //  %% load data from file to variable data, data=(Y,lambda)
            //  i=find(data(:,1)>0.15)
            //  data = data = data(i,:)
            //  c = polyfit(data(:,1),data(:,2),3)
            //  %% To plot
            //  yy = 0:0.01:1
            //  fit = c(1)*yy.^3 + c(2)*yy.^2 + c(3)*yy + c(4);
            //  plot(data(:,1),data(:,2),'*')
            //  hold on
            //  plot(yy,fit,'r')
            double lam_fit = -0.026920*Y*Y*Y + 0.066128*Y*Y - 0.056929*Y + 0.978174;

            if(Y>0.6)
            {
                double error = 0.0004;
                TS_ASSERT_DELTA(lambda[i], lam_fit, error);
            }
            else if(Y>0.15)
            {
                double error = 0.0015;
                TS_ASSERT_DELTA(lambda[i], lam_fit, error);
            }

            //// don't delete:
            //std::cout << quad_points.Get(i)(0) << " " << quad_points.Get(i)(1) << " " << lambda[i] << "\n";
        }

        // hardcoded test
        TS_ASSERT_DELTA(lambda[34], 0.9753, 1e-4);
    }


    // cover all other contraction model options which are allowed but not been used in a test 
    // so far (or in TestExplicitCardiacMechanicsAssembler)
    void TestCoverage() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1, 1);

        MooneyRivlinMaterialLaw<2> law(1);
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        ImplicitCardiacMechanicsAssembler<2> expl_solver1(KERCHOFFS2003,&mesh,"",fixed_nodes,&law);
        ImplicitCardiacMechanicsAssembler<2> expl_solver2(NONPHYSIOL3,&mesh,"",fixed_nodes,&law);

        // call with TS_ASSERT_THROWS_CONTAINS with any disallowed contraction models here:
    }
};

#endif /*TESTIMPLICITCARDIACMECHANICSASSEMBLER_HPP_*/
