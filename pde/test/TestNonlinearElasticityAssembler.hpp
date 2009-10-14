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


#ifndef TESTNONLINEARELASTICITYASSEMBLER_HPP_
#define TESTNONLINEARELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "NonlinearElasticityAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "MooneyRivlinMaterialLaw.hpp"

double MATERIAL_PARAM = 0.05;
double ALPHA = 0.2;

// Body force corresponding to the deformation
// x = X+0.5*alpha*X^2, y=Y/(1+alpha*X), with p=2c
c_vector<double,2> MyBodyForce(c_vector<double,2>& X)
{
    assert(X(0)>=0 && X(0)<=1 && X(1)>=0 && X(1)<=1);

    c_vector<double,2> body_force;
    double lam = 1+ALPHA*X(0);
    body_force(0) = -2*MATERIAL_PARAM * ALPHA;
    body_force(1) = -2*MATERIAL_PARAM * 2*ALPHA*ALPHA*X(1)/(lam*lam*lam);
    return body_force;
}

// Surface traction on three sides of a cube, corresponding to
// x = X+0.5*alpha*X^2, y=Y/(1+alpha*X), with p=2c
c_vector<double,2> MyTraction(c_vector<double,2>& X)
{
    c_vector<double,2> traction = zero_vector<double>(2);

    double lam = 1+ALPHA*X(0);
    if (X(0)==1)
    {
        traction(0) =  2*MATERIAL_PARAM * (lam - 1.0/lam);
        traction(1) = -2*MATERIAL_PARAM * X(1)*ALPHA/(lam*lam);
    }
    else if (X(1)==0)
    {
        traction(0) =  2*MATERIAL_PARAM * X(1)*ALPHA/(lam*lam);
        traction(1) = -2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else if (X(1)==1)
    {
        traction(0) = -2*MATERIAL_PARAM * X(1)*ALPHA/(lam*lam);
        traction(1) =  2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else
    {
        NEVER_REACHED;
    }
    return traction;
}



class TestNonlinearElasticityAssembler : public CxxTest::TestSuite
{
public:
    void TestAssembleSystem() throw (Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);
        ExponentialMaterialLaw<2> law(2,3);
        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        NonlinearElasticityAssembler<2> assembler(&mesh,
                                                  &law,
                                                  zero_vector<double>(2),
                                                  1.0,
                                                  "",
                                                  fixed_nodes);
        assembler.AssembleSystem(true, true);

        ///////////////////////////////////////////////////////////////////
        // test whether residual vector is currently zero (as
        // current solution should have been initialised to u=0, p=p0
        ///////////////////////////////////////////////////////////////////
        ReplicatableVector rhs_vec(assembler.mpLinearSystem->rGetRhsVector());
        TS_ASSERT_EQUALS( rhs_vec.GetSize(), 2U*289U+81U );
        for (unsigned i=0; i<rhs_vec.GetSize(); i++)
        {
            TS_ASSERT_DELTA(rhs_vec[i], 0.0, 1e-12);
        }

        ///////////////////////////////////////////////////////////////////
        // compute numerical jacobian and compare with analytic jacobian
        // (about u=0, p=p0)
        ///////////////////////////////////////////////////////////////////
        unsigned num_dofs = rhs_vec.GetSize();
        double h = 1e-6;

        int lo, hi;
        MatGetOwnershipRange(assembler.mpLinearSystem->rGetLhsMatrix(), &lo, &hi);

        for (unsigned j=0; j<num_dofs; j++)
        {
            assembler.mCurrentSolution.clear();
            assembler.FormInitialGuess();
            assembler.mCurrentSolution[j] += h;

            assembler.AssembleSystem(true, false);

            ReplicatableVector perturbed_rhs( assembler.mpLinearSystem->rGetRhsVector() );

            for (unsigned i=0; i<num_dofs; i++)
            {
                if ((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = assembler.mpLinearSystem->GetMatrixElement(i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec[i])/h;
                    if ((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-2);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-2);
                    }
                }
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);

        //////////////////////////////////////////////////////////
        // compare numerical and analytic jacobians again, this
        // time using a non-zero displacement, u=lambda x, v = mu y
        // (lambda not equal to 1/nu), p = p0.
        //////////////////////////////////////////////////////////
        double lambda = 1.2;
        double mu = 1.0/1.3;

        assembler.mCurrentSolution.clear();
        assembler.FormInitialGuess();
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            assembler.mCurrentSolution[2*i]   = (lambda-1)*mesh.GetNode(i)->rGetLocation()[0];
            assembler.mCurrentSolution[2*i+1] = (mu-1)*mesh.GetNode(i)->rGetLocation()[1];
        }

        assembler.AssembleSystem(true, true);
        ReplicatableVector rhs_vec2(assembler.mpLinearSystem->rGetRhsVector());

        h=1e-8; // needs to be smaller for this one

        for (unsigned j=0; j<num_dofs; j++)
        {
            assembler.mCurrentSolution[j] += h;
            assembler.AssembleSystem(true, false);

            ReplicatableVector perturbed_rhs( assembler.mpLinearSystem->rGetRhsVector() );

            for (unsigned i=0; i<num_dofs; i++)
            {
                if ((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = assembler.mpLinearSystem->GetMatrixElement(i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec2[i])/h;
                    if ((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-2);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-2);
                    }
                }
            }
            assembler.mCurrentSolution[j] -= h;
        }
    }


    // A test where the solution should be zero displacement
    // It mainly tests that the initial guess was set up correctly to
    // the final correct solution, ie u=0, p=zero_strain_pressure (!=0)
    void TestWithZeroDisplacement() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        double c1 = 3.0;
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(c1);

        std::vector<unsigned> fixed_nodes;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if ( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }

        NonlinearElasticityAssembler<2> assembler(&mesh,
                                                  &mooney_rivlin_law,
                                                  zero_vector<double>(2),
                                                  1.0,
                                                  "",
                                                  fixed_nodes);

        // for coverage
        TS_ASSERT_THROWS_THIS(assembler.SetWriteOutput(true),
                "Can\'t write output if no output directory was given in constructor");
        assembler.SetWriteOutput(false);

        assembler.Solve();
        TS_ASSERT_EQUALS(assembler.GetNumNewtonIterations(), 0u);

        // get deformed position
        std::vector<c_vector<double,2> >& r_deformed_position
            = assembler.rGetDeformedPosition();

        for (unsigned i=0; i<r_deformed_position.size(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], r_deformed_position[i](0), 1e-8);
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[1], r_deformed_position[i](1), 1e-8);
        }

        // check the final pressure
        std::vector<double>& r_pressures = assembler.rGetPressures();
        TS_ASSERT_EQUALS(r_pressures.size(), mesh.GetNumVertices());
        for (unsigned i=0; i<r_pressures.size(); i++)
        {
            TS_ASSERT_DELTA(r_pressures[i], 2*c1, 1e-6);
        }
    }

    void TestSettingUpHeterogeneousProblem() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // two element quad mesh on the square
        QuadraticMesh<2> mesh(1.0, 1.0, 1, 1);

        MooneyRivlinMaterialLaw<2> law_1(1.0);
        MooneyRivlinMaterialLaw<2> law_2(5.0);
        std::vector<AbstractIncompressibleMaterialLaw<2>*> laws;
        laws.push_back(&law_1);
        laws.push_back(&law_2);

        std::vector<unsigned> fixed_nodes;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if ( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }

        NonlinearElasticityAssembler<2> assembler(&mesh,
                                                  laws,
                                                  zero_vector<double>(2),
                                                  1.0,
                                                  "",
                                                  fixed_nodes);

        TS_ASSERT_EQUALS(assembler.mMaterialLaws.size(), 2u);
        TS_ASSERT_DELTA(assembler.mMaterialLaws[0]->GetZeroStrainPressure(), 2.0, 1e-6);
        TS_ASSERT_DELTA(assembler.mMaterialLaws[1]->GetZeroStrainPressure(), 10.0, 1e-6);

        unsigned num_nodes = 9;
        // pressure for node 0 (in elem 0)
        TS_ASSERT_DELTA(assembler.mCurrentSolution[2*num_nodes + 0], 2.0, 1e-6);
        // pressure for node 3 (in elem 1)
        TS_ASSERT_DELTA(assembler.mCurrentSolution[2*num_nodes + 3], 10.0, 1e-6);
    }

    void TestSolve() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        MooneyRivlinMaterialLaw<2> law(0.02);
        c_vector<double,2> body_force;
        body_force(0) = 0.06;
        body_force(1) = 0.0;

        std::vector<unsigned> fixed_nodes;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if ( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }

        NonlinearElasticityAssembler<2> assembler(&mesh,
                                                  &law,
                                                  body_force,
                                                  1.0,
                                                  "simple_nonlin_elas",
                                                  fixed_nodes);

        assembler.Solve();

        std::vector<c_vector<double,2> >& r_solution = assembler.rGetDeformedPosition();

        double xend = 1.17199;
        double yend = 0.01001;

        ///////////////////////////////////////////////////////////
        // compare the solution at the corners with the values
        // obtained using the dealii finite elasticity assembler
        //
        // Results have been visually checked to see they agree
        // (they do, virtually or completely overlapping.
        ///////////////////////////////////////////////////////////

        // bottom lhs corner should still be at (0,0)
        assert( fabs(mesh.GetNode(0)->rGetLocation()[0] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(0)->rGetLocation()[1] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[0](0), 0.0, 1e-6 );
        TS_ASSERT_DELTA( r_solution[0](1), 0.0, 1e-6 );

        // top lhs corner should still be at (0,1)
        assert( fabs(mesh.GetNode(3)->rGetLocation()[0] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(3)->rGetLocation()[1] - 1) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[3](0), 0.0, 1e-6 );
        TS_ASSERT_DELTA( r_solution[3](1), 1.0, 1e-6 );

        // DEALII value for bottom rhs corner is (1.17199,0.01001)
        assert( fabs(mesh.GetNode(1)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(1)->rGetLocation()[1] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[1](0), xend, 1e-3 );
        TS_ASSERT_DELTA( r_solution[1](1), yend, 1e-3 );

        // DEALII value for top rhs corner is (1.17199,0.98999)
        assert( fabs(mesh.GetNode(2)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(2)->rGetLocation()[1] - 1) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[2](0), xend,   1e-3 );
        TS_ASSERT_DELTA( r_solution[2](1), 1-yend, 1e-3 );
    }



    /**
     *  Solve a problem with non-zero dirichlet boundary conditions
     *  and non-zero tractions. THIS TEST COMPARES AGAINST AN EXACT SOLUTION.
     *
     *  Choosing the deformation x=X/lambda, y=lambda*Y, with a
     *  Mooney-Rivlin material, then
     *   F = [1/lam 0; 0 lam], T = [2*c1-p*lam^2, 0; 0, 2*c1-p/lam^2],
     *   sigma = [2*c1/lam^2-p, 0; 0, 2*c1*lam^2-p].
     *  Choosing p=2*c1*lam^2, then sigma = [2*c1/lam^2-p 0; 0 0].
     *  The surface tractions are then
     *   TOP and BOTTOM SURFACE: 0
     *   RHS: s = SN = J*invF*sigma*N = [lam 0; 0 1/lam]*sigma*[1,0]
     *          = [2*c1(1/lam-lam^3), 0]
     *
     *  So, we have to specify displacement boundary conditions (y=lam*Y) on
     *  the LHS (X=0), and traction bcs (s=the above) on the RHS (X=1), and can
     *  compare the computed displacement and pressure against the true solution.
     *
     */
    void TestSolveWithNonZeroBoundaryConditions() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        double lambda = 0.85;
        double c1 = 0.02;
        c_vector<double,2> body_force = zero_vector<double>(2);
        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0, 1.0, num_elem, num_elem);
        MooneyRivlinMaterialLaw<2> law(c1);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > locations;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if ( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
                c_vector<double,2> new_position;
                new_position(0) = 0;
                new_position(1) = lambda*mesh.GetNode(i)->rGetLocation()[1];
                locations.push_back(new_position);
            }
        }

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        c_vector<double,2> traction;
        traction(0) = 2*c1*(pow(lambda,-1) - lambda*lambda*lambda);
        traction(1) = 0;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[0] - 1.0)<1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }
        assert(boundary_elems.size()==num_elem);

        NonlinearElasticityAssembler<2> assembler(&mesh,
                                                  &law,
                                                  body_force,
                                                  1.0,
                                                  "nonlin_elas_non_zero_bcs",
                                                  fixed_nodes,
                                                  &locations);

        assembler.SetSurfaceTractionBoundaryConditions(boundary_elems, tractions);

        assembler.Solve();

        std::vector<c_vector<double,2> >& r_solution = assembler.rGetDeformedPosition();

        for (unsigned i=0; i<fixed_nodes.size(); i++)
        {
            unsigned index = fixed_nodes[i];
            TS_ASSERT_DELTA(r_solution[index](0), locations[i](0), 1e-8);
            TS_ASSERT_DELTA(r_solution[index](1), locations[i](1), 1e-8);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double exact_x = (1.0/lambda)*mesh.GetNode(i)->rGetLocation()[0];
            double exact_y = lambda*mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-5 );
            TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-5 );
        }

        for (unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
            TS_ASSERT_DELTA( assembler.rGetPressures()[i], 2*c1*lambda*lambda, 1e-6 );
        }
    }

    /**
     *  Test with functional (rather than constant) body force and surface traction, against a known
     *  solution. Since a non-zero body force is used here and a known solution, this is the MOST
     *  IMPORTANT TEST.
     *
     *  Choose x=X+0.5*alpha*X*X, y=Y/(1+alpha*X), p=2c, then F has determinant 1, and S can be shown to
     *  be, where lam = 1 +alpha X (ie dx/dX)
     *    S = 2c[lam-1/lam,   -Y*alpha*lam^{-2}; -Y*alpha*lam^{-2}, 1/lam - lam]
     *  in which case the required body force and surface traction can be computed to be
     *
     *  b = 2c/density [ -alpha, -2*Y*alpha^2 * lam^{-3} ]
     *  s = 2c[lam-1/lam, -Y*alpha/(lam^2)] on X=1
     *  s = 2c[0, lam - 1/lam]              on Y=0
     *  s = 2c[Y*alpha/lam^2, 1/lam - lam]  on Y=1
     *
     */
    void TestWithFunctionalData() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools
        MechanicsEventHandler::Reset();

        c_vector<double,2> body_force = zero_vector<double>(2);

        unsigned num_elem = 5;
        QuadraticMesh<2> mesh(1.0, 1.0, num_elem, num_elem);

        MooneyRivlinMaterialLaw<2> law(MATERIAL_PARAM);

        std::vector<unsigned> fixed_nodes;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if ( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            // get all boundary elems except those on X=0
            if (fabs((*iter)->CalculateCentroid()[0])>1e-6)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
            }
        }
        assert(boundary_elems.size()==3*num_elem);

        NonlinearElasticityAssembler<2> assembler(&mesh,
                                                  &law,
                                                  body_force,
                                                  1.0,
                                                  "nonlin_elas_functional_data",
                                                  fixed_nodes);

        assembler.SetFunctionalBodyForce(MyBodyForce);
        assembler.SetFunctionalTractionBoundaryCondition(boundary_elems, MyTraction);

        assembler.Solve();

        // matrix might have (small) errors introduced if this fails
        TS_ASSERT_EQUALS(assembler.GetNumNewtonIterations(), 3u);

        std::vector<c_vector<double,2> >& r_solution = assembler.rGetDeformedPosition();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double X = mesh.GetNode(i)->rGetLocation()[0];
            double Y = mesh.GetNode(i)->rGetLocation()[1];

            double exact_x = X + 0.5*ALPHA*X*X;
            double exact_y = Y/(1+ALPHA*X);

            TS_ASSERT_DELTA(r_solution[i](0), exact_x, 1e-4);
            TS_ASSERT_DELTA(r_solution[i](1), exact_y, 1e-4);
        }

        for (unsigned i=0; i<assembler.rGetPressures().size(); i++)
        {
            TS_ASSERT_DELTA( assembler.rGetPressures()[i]/(2*MATERIAL_PARAM), 1.0, 1e-3);
        }

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }
};

#endif /*TESTNONLINEARELASTICITYASSEMBLER_HPP_*/
