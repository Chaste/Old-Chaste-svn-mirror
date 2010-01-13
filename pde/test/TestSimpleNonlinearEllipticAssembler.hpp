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
#ifndef _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_



#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>
//#include <iostream>

#include "SimpleNonlinearEllipticAssembler.hpp"
#include "SimplePetscNonlinearSolver.hpp"

#include "BoundaryConditionsContainer.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "ConstBoundaryCondition.hpp"

#include "NonlinearEquationPde.hpp"
#include "NonlinearEquation2Pde.hpp"
#include "NonlinearEquation3Pde.hpp"
#include "NonlinearEquation4Pde.hpp"
#include "NonlinearEquation5Pde.hpp"
#include "Example2DNonlinearEllipticPde.hpp"
#include "NonlinearLinearEquation.hpp"
#include "ExampleNasty2dNonlinearEllipticPde.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"

/**
 * For use in TestSimpleNonlinearEllipticAssembler::Test2dOnUnitSquare.
 */
double bc_x1_func(const ChastePoint<2>& p)
{
    return 2*(2+p[1]*p[1]);
}
/**
 * For use in TestSimpleNonlinearEllipticAssembler::Test2dOnUnitSquare.
 */
double bc_y1_func(const ChastePoint<2>& p)
{
    return 2*(2+p[0]*p[0]);
}

/**
 * For use in TestSimpleNonlinearEllipticAssembler::TestNasty2dEquationOnUnitSquare.
 */
double bc_x1_func2(const ChastePoint<2>& p)
{
    return sin(2)*(sin(1)*sin(1)+1+p[1]*p[1]);
}
/**
 * For use in TestSimpleNonlinearEllipticAssembler::TestNasty2dEquationOnUnitSquare.
 */
double bc_y1_func2(const ChastePoint<2>& p)
{
    return 2*(2+sin(p[0])*sin(p[0]));
}

/**
 * For use in TestSimpleNonlinearEllipticAssembler::TestWithHeatEquation2DAndNeumannBCs
 */
double one_bc(const ChastePoint<2>& p)
{
    return p[1];
}


class TestSimpleNonlinearEllipticAssembler : public CxxTest::TestSuite
{
public:
    void TestAssembleResidual()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/practical1_1d_mesh");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        //Adding Dirichlet BC at node 0
        double DirichletBCValue = 5.0;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(DirichletBCValue);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), pBoundaryCondition);

        // adding von Neumann BC at the last node
        double VonNeumannBCValue = 9.0;
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(VonNeumannBCValue);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--; // to be consistent with c++ :))), GetBoundaryElementIteratorEnd points to one element passed it
        bcc.AddNeumannBoundaryCondition(*iter,pBoundaryCondition1);

        // initialize 'solution' vector
        double initial_guess_value = 1.0;
        double h = 0.01;
        Vec solution = PetscTools::CreateVec(mesh.GetNumNodes(), initial_guess_value);

        NonlinearEquationPde<1> pde;
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);

        Vec residual;
        VecDuplicate(solution, &residual);

        assembler.PrepareForSolve();
        assembler.AssembleResidual(solution, residual);

        ReplicatableVector residual_repl(residual);
        TS_ASSERT(fabs(residual_repl[0] + DirichletBCValue - initial_guess_value) < 0.001);
        TS_ASSERT(fabs(residual_repl[1] + h) < 0.001);
        TS_ASSERT(fabs(residual_repl[mesh.GetNumNodes()-1] + VonNeumannBCValue + h/2) < 0.001);

        VecDestroy(residual);
        VecDestroy(solution);
    }

    void TestNumericalAgainstAnalyticJacobian()
    {
        PetscInt n = 11;  // Mesh size
        Mat numerical_jacobian;
        PetscTools::SetupMat(numerical_jacobian, n, n, (MatType) MATSEQDENSE);

        Mat analytic_jacobian;
        PetscTools::SetupMat(analytic_jacobian,  n, n, (MatType) MATSEQDENSE);

        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);


        // assembler
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);
        assembler.PrepareForSolve();

        // cover VerifyJacobian
        TS_ASSERT( assembler.VerifyJacobian(1e-3) );

        // Set up initial solution guess for residuals
        std::vector<double> init_guess(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            init_guess[i] = -0.01*i*i;
        }
        Vec initial_guess = PetscTools::CreateVec(init_guess);


        int errcode = assembler.AssembleJacobianNumerically(initial_guess, &numerical_jacobian);
        TS_ASSERT_EQUALS(errcode, 0);

        assembler.mUseAnalyticalJacobian = true; // can access the member variable as this class is a friend
        errcode = assembler.AssembleJacobian(initial_guess, &analytic_jacobian);

        TS_ASSERT_EQUALS(errcode, 0);
        MatAssemblyBegin(numerical_jacobian, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(numerical_jacobian, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(analytic_jacobian, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(analytic_jacobian, MAT_FINAL_ASSEMBLY);

        //TS_TRACE("Numerical:");
        //MatView(numerical_jacobian, 0);
        //TS_TRACE("Analytical:");
        //MatView(analytic_jacobian, 0);

        PetscScalar numerical_array[n*n], analytic_array[n*n];
        PetscInt row_ids[n], col_ids[n];
        int lo, hi;
        MatGetOwnershipRange(numerical_jacobian, &lo, &hi);
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            row_ids[local_index] = global_index;
        }
        for (int i=0; i<n; i++)
        {
            col_ids[i] = i;
        }

        // Check matrices are the same, to within numerical error tolerance
        MatGetValues(numerical_jacobian, hi-lo, row_ids, n, col_ids, numerical_array);
        MatGetValues(analytic_jacobian, hi-lo, row_ids, n, col_ids, analytic_array);
        for (int local_index=0; local_index<hi-lo; local_index++)
        {
            for (int j=0; j<n; j++)
            {
                TS_ASSERT_DELTA(numerical_array[local_index*n+j],
                                analytic_array[local_index*n+j], 0.001);
            }
        }
        VecDestroy(initial_guess);
        MatDestroy(numerical_jacobian);
        MatDestroy(analytic_jacobian);
    }

    void TestWithHeatEquation1D()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);

        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);

        // Set up initial guess
        Vec initial_guess = PetscTools::CreateVec(mesh.GetNumNodes(),1.0);

        Vec answer = assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = sqrt(x*(1-x));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.001);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void TestWithHeatEquation1DAndNeumannBCs()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        // u(0) = 0
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        // u(1)*u'(1) = 1
        p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        // Nonlinear assembler to use
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = assembler.CreateConstantInitialGuess(0.25);

        // Set no. of gauss points to use
        assembler.SetNumberOfQuadraturePointsPerDimension(3);

        // Solve the PDE
        Vec answer = assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = sqrt(x*(4-x));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.001);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void TestWithHeatEquation1D2()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquation2Pde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        p_boundary_condition = new ConstBoundaryCondition<1>(exp(1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);

        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);

        // Set up initial Guess
        std::vector<double> init_guess(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            init_guess[i] = 1.0+0.01*i*i;
        }
        Vec initial_guess = PetscTools::CreateVec(init_guess);

        Vec answer = assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = exp(0.5*(3.0*x-x*x));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.001);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void TestWithHeatEquation1D3()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquation3Pde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(sqrt(2.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);

        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc,3);

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateVec(mesh.GetNumNodes());
        for (unsigned global_index=0; global_index<mesh.GetNumNodes(); global_index++)
        {
            VecSetValue(initial_guess, global_index, (1.5-0.15*global_index), INSERT_VALUES);
        }
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);

        Vec answer = assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = sqrt(2.0*(exp(-x)-x*exp(-1.0)));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.001);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void TestWithHeatEquation1D4()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquation4Pde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        // u(1) = exp(1.0)
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(exp(-1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);
        // u(0)^2*u'(0) = 0.0
        p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);

        // Set up initial Guess
        std::vector<double> init_guess(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x1=0.1*(double)(i);
            init_guess[i] =  0.35*(1-x1*x1);
        }
        Vec initial_guess = PetscTools::CreateVec(init_guess);

        Vec answer = assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = x*exp(-x);
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void TestWithHeatEquation1D5()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquation5Pde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        // u(1) = exp(-1.0)
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(exp(-1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);
        // u(0)^2*u'(0) = -1.0
        // Note that we specify 1 as the value, since figuring out which direction
        // the normal is in is hard in 1D.
        p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateVec(mesh.GetNumNodes());
        double x1;
        for (unsigned global_index=0; global_index<mesh.GetNumNodes(); global_index++)
        {
            x1=0.1*(double)(global_index);
            VecSetValue(initial_guess, global_index, 0.35*(1-x1*x1), INSERT_VALUES);
        }
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);

        Vec answer = assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = exp(-x);
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void TestWithHeatEquation1DAndNeumannBCs2()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        // u(1) = sqrt(3)
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(sqrt(3));
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(10), p_boundary_condition);
        // u(0)*u'(0) = 2
        // Note that we specify -2 as the value, since figuring out which direction
        // the normal is in is hard in 1D.
        p_boundary_condition = new ConstBoundaryCondition<1>(-2.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);

        // cover the bad size exception
        Vec badly_sized_init_guess = PetscTools::CreateVec(1,1.0); // size=1
        TS_ASSERT_THROWS_THIS( assembler.Solve(badly_sized_init_guess, true),
                "Size of initial guess vector, 1, does not match size of problem, 11" );

        // Set up initial Guess
        Vec initial_guess = assembler.CreateConstantInitialGuess(1.0);

        // This problem seems unusally sensitive to the initial guess. Various other
        // choices failed to converge.
        Vec answer = assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = sqrt(x*(4-x));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.001);
        }

        VecDestroy(badly_sized_init_guess);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void TestHeatEquationWithNeumannOnUnitDisc()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearLinearEquation<2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        // du/dn = -0.5 on r=1
        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_boundary_condition;
        p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition);

        SimpleNonlinearEllipticAssembler<2,2> assembler(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = assembler.CreateConstantInitialGuess(1.0);

        Vec answer = assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            c_vector<double, 2> r;
            r(0) = mesh.GetNode(i)->GetPoint()[0];
            r(1) = mesh.GetNode(i)->GetPoint()[1];
            double u = -0.25 * inner_prod(r, r) + 2.25;
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void TestWithHeatEquation2DAndNeumannBCs()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearEquationPde<2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        // u(y=0) = 0
        ConstBoundaryCondition<2>* zero_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
        while (node_iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            if (fabs((*node_iter)->GetPoint()[1]) < 1e-12)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, zero_boundary_condition);
            }
            node_iter++;
        }

        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        FunctionalBoundaryCondition<2>* one_boundary_condition = new FunctionalBoundaryCondition<2>(&one_bc);
        AbstractBoundaryCondition<2>* p_boundary_condition;
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            double y = (*iter)->GetNodeLocation(0,1);
            if (fabs(y-1.0) < 1e-12)
            {
                // u(y=1)*u'(y=1) = 1
                p_boundary_condition = one_boundary_condition;
            }
            else
            {
                // No flux across left & right
                p_boundary_condition = zero_boundary_condition;
            }

            bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

            iter++;
        }

        SimpleNonlinearEllipticAssembler<2,2> assembler(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = assembler.CreateConstantInitialGuess(0.25);

        // solve
        Vec answer = assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = sqrt(y*(4-y));
            TS_ASSERT_DELTA(answer_repl[i], u, 0.15);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void Test2dOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        Example2DNonlinearEllipticPde pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition;
        TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
        while (node_iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*node_iter)->GetPoint()[0];
            double y = (*node_iter)->GetPoint()[1];
            p_boundary_condition = NULL;
            if (fabs(x) < 1e-12)
            {
                // On x=0, u=1+y^2
                p_boundary_condition = new ConstBoundaryCondition<2>(1 + y*y);
            }
            else if (fabs(y) < 1e-12)
            {
                // On y=0, u=1+x^2
                p_boundary_condition = new ConstBoundaryCondition<2>(1 + x*x);
            }
            if (p_boundary_condition)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition);
            }

            node_iter++;
        }
        FunctionalBoundaryCondition<2>* p_functional_bc;
        TetrahedralMesh<2,2>::BoundaryElementIterator elt_iter = mesh.GetBoundaryElementIteratorBegin();
        while (elt_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            double x = (*elt_iter)->GetNodeLocation(0,0);
            double y = (*elt_iter)->GetNodeLocation(0,1);
            p_functional_bc = NULL;
            if (fabs(y-1.0) < 1e-12)
            {
                // On y=1, Dgradu_dot_n = 2(2+x^2)
                p_functional_bc = new FunctionalBoundaryCondition<2>(&bc_y1_func);
            }
            else if (fabs(x-1.0) < 1e-12)
            {
                // On x=1, Dgradu_dot_n = 2(2+y^2)
                p_functional_bc = new FunctionalBoundaryCondition<2>(&bc_x1_func);
            }
            if (p_functional_bc)
            {
                bcc.AddNeumannBoundaryCondition(*elt_iter, p_functional_bc);
            }

            elt_iter++;
        }

        SimpleNonlinearEllipticAssembler<2,2> assembler(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = assembler.CreateConstantInitialGuess(4.0);


        // Numerical Jacobian
        Vec answer = assembler.Solve(initial_guess, false);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = 1 + x*x + y*y;
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        VecDestroy(answer);

        // Analytical Jacobian
        answer=assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl2(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = 1 + x*x + y*y;
            TS_ASSERT_DELTA(answer_repl2[i], u, 0.01);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }

    void TestNasty2dEquationOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        ExampleNasty2dNonlinearEllipticPde pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition;
        TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
        while (node_iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*node_iter)->GetPoint()[0];
            double y = (*node_iter)->GetPoint()[1];
            p_boundary_condition = NULL;
            if (fabs(x) < 1e-12)
            {
                // On x=0, u=1+y^2
                p_boundary_condition = new ConstBoundaryCondition<2>(1 + y*y);
            }
            else if (fabs(y) < 1e-12)
            {
                // On y=0, u=1+sin^2(x)
                p_boundary_condition = new ConstBoundaryCondition<2>(1 + sin(x)*sin(x));
            }
            if (p_boundary_condition)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition);
            }

            node_iter++;
        }
        FunctionalBoundaryCondition<2>* p_functional_bc;
        TetrahedralMesh<2,2>::BoundaryElementIterator elt_iter = mesh.GetBoundaryElementIteratorBegin();
        while (elt_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            double x = (*elt_iter)->GetNodeLocation(0,0);
            double y = (*elt_iter)->GetNodeLocation(0,1);
            p_functional_bc = NULL;
            if (fabs(y-1.0) < 1e-12)
            {
                // On y=1, Dgradu_dot_n = 2(2+sin^2(x))
                p_functional_bc = new FunctionalBoundaryCondition<2>(&bc_y1_func2);
            }
            else if (fabs(x-1.0) < 1e-12)
            {
                // On x=1, Dgradu_dot_n = sin(2)(sin^2(1)+1+y^2)
                p_functional_bc = new FunctionalBoundaryCondition<2>(&bc_x1_func2);
            }
            if (p_functional_bc)
            {
                bcc.AddNeumannBoundaryCondition(*elt_iter, p_functional_bc);
            }

            elt_iter++;
        }

        SimpleNonlinearEllipticAssembler<2,2> assembler(&mesh, &pde, &bcc);

        // Set up initial Guess
        Vec initial_guess = assembler.CreateConstantInitialGuess(4.0);

        // Numerical Jacobian
        Vec answer = assembler.Solve(initial_guess, false);
        ReplicatableVector answer_repl(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = 1 + sin(x)*sin(x) + y*y;
            TS_ASSERT_DELTA(answer_repl[i], u, 0.01);
        }

        VecDestroy(answer);

        // Analytical Jacobian
        answer=assembler.Solve(initial_guess, true);
        ReplicatableVector answer_repl2(answer);

        // Check result
        for (unsigned i=0; i<answer_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = 1 + sin(x)*sin(x) + y*y;
            TS_ASSERT_DELTA(answer_repl2[i], u, 0.01);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
};

#endif //_TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
