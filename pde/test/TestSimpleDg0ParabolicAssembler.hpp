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
#ifndef _TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_

/**
 * TestSimpleDg0ParabolicAssembler.hpp
 *
 * Test suite for the Dg0ParabolicAssembler class.
 *
 * Tests the class for the solution of parabolic pdes in 1D, 2D and 3D with and
 * without source terms with neumann and dirichlet booundary conditions.
 */

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "ParallelColumnDataWriter.hpp"
#include "TrianglesMeshReader.hpp"
#include "FemlabMeshReader.hpp"
#include "HeatEquation.hpp"
#include "HeatEquationWithSourceTerm.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"


class TestSimpleDg0ParabolicAssembler : public CxxTest::TestSuite
{
public:

    void TestExceptionalBehaviour()
    {
        // Assembler
        SimpleDg0ParabolicAssembler<1,1, true> assembler(NULL,NULL,NULL);

        // start > end
        TS_ASSERT_THROWS_ANYTHING(assembler.SetTimes(1.0, 0.0, 0.01));

        // dt = 0
        TS_ASSERT_THROWS_ANYTHING(assembler.SetTimes(0.0, 1.0, 0.0));
    }

    /// test 1D problem
    void TestSimpleDg0ParabolicAssembler1DZeroDirich()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<1> pde;

        // Boundary conditions - zero dirichlet at first and last node;
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1> *p_boundary_condition =
            new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode( mesh.GetNumNodes()-1 ), p_boundary_condition);

        // Assembler
        SimpleDg0ParabolicAssembler<1,1,true> assembler(&mesh,&pde,&bcc);

        // Initial condition, u(0,x) = sin(x*pi);
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = sin(x*M_PI);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);

        double t_end = 0.1;
        assembler.SetTimes(0, t_end, 0.01);
        assembler.SetInitialCondition(initial_condition);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = e^{-t*pi*pi} sin(x*pi), t=1
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = exp(-0.1*M_PI*M_PI)*sin(x*M_PI);
            TS_ASSERT_DELTA(result_repl[i], u, 0.1);
        }

        VecDestroy(initial_condition);
        VecDestroy(result);
    }


    void TestSimpleDg0ParabolicAssembler1DZeroDirichWithSourceTerm()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquationWithSourceTerm<1> pde;

        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1> *p_boundary_condition =
            new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        p_boundary_condition = new ConstBoundaryCondition<1>(-0.5);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode( mesh.GetNumNodes()-1 ), p_boundary_condition);

        // Assembler
        SimpleDg0ParabolicAssembler<1,1,true> assembler(&mesh,&pde,&bcc);

        // initial condition, u(0,x) = sin(x*pi)+0.5*x*x;
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = sin(x*M_PI)-0.5*x*x;
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);

        double t_end = 0.1;
        assembler.SetTimes(0, t_end, 0.01);
        assembler.SetInitialCondition(initial_condition);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = e^{-t*pi*pi} sin(x*pi) + 0.5*x^2, t=1
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = exp(-0.1*M_PI*M_PI)*sin(x*M_PI)-0.5*x*x;
            TS_ASSERT_DELTA(result_repl[i], u, 0.1);
        }

        VecDestroy(initial_condition);
        VecDestroy(result);
    }

    void TestSimpleDg0ParabolicAssemblerNonzeroNeumannCondition()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<1> pde;

        // Boundary conditions  u(0)=0, u'(1)=1
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1> *p_boundary_condition = new ConstBoundaryCondition<1>(0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);

        ConstBoundaryCondition<1> *p_neumann_boundary_condition =
            new ConstBoundaryCondition<1>(1.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);

        // Assembler
        SimpleDg0ParabolicAssembler<1,1,true> assembler(&mesh,&pde,&bcc);

        // initial condition;
        const double PI_over_2 = M_PI/2.0;
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = x + sin(PI_over_2 * x);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);

        // set time and initial condition
        assembler.SetTimes(0, 0.5, 0.01);
        assembler.SetInitialCondition(initial_condition);

        // solve
        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // check result
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = x + exp(-0.5*PI_over_2*PI_over_2)*sin(x*PI_over_2);
            TS_ASSERT_DELTA(result_repl[i], u, 0.01);
        }

        VecDestroy(initial_condition);
        VecDestroy(result);
    }


    void TestSimpleDg0ParabolicAssembler2DZeroDirich()
    {
        // read mesh on [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        // Boundary conditions - zero dirichlet everywhere on boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Assembler
        SimpleDg0ParabolicAssembler<2,2,true> assembler(&mesh,&pde,&bcc);

        // initial condition;
        // choose initial condition sin(x*pi)*sin(y*pi) as this is an eigenfunction of
        // the heat equation.
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(x*M_PI)*sin(y*M_PI);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);

        // Solve
        double t_end = 0.1;
        assembler.SetTimes(0, t_end, 0.001);
        assembler.SetInitialCondition(initial_condition);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // check result
        // Solution should be u = e^{-2*t*pi*pi} sin(x*pi)*sin(y*pi), t=1
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = exp(-2*t_end*M_PI*M_PI)*sin(x*M_PI)*sin(y*M_PI);
            TS_ASSERT_DELTA(result_repl[i], u, 0.01);
        }

        VecDestroy(initial_condition);
        VecDestroy(result);
    }


    // test 2D problem
    void TestSimpleDg0ParabolicAssembler2DZeroDirichWithSourceTerm()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquationWithSourceTerm<2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();

        while (iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            ConstBoundaryCondition<2> *p_dirichlet_boundary_condition =
                new ConstBoundaryCondition<2>(-0.25*(x*x+y*y));
            bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            iter++;
        }

        // Assembler
        SimpleDg0ParabolicAssembler<2,2,true> assembler(&mesh,&pde,&bcc);

        // initial condition, u(0,x) = sin(x*pi)*sin(y*pi)-0.25*(x^2+y^2);
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(x*M_PI)*sin(y*M_PI)-0.25*(x*x+y*y);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);

        double t_end = 0.1;
        assembler.SetTimes(0, t_end, 0.001);
        assembler.SetInitialCondition(initial_condition);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // check result
        // Solution should be u = e^{-t*2*pi*pi} sin(x*pi) sin(y*pi) - 0.25(x^2+y^2), t=0.1
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = exp(-0.1*2*M_PI*M_PI)*sin(x*M_PI)*sin(y*M_PI)-0.25*(x*x+y*y);
            TS_ASSERT_DELTA(result_repl[i], u, 0.05);
        }

        VecDestroy(initial_condition);
        VecDestroy(result);
    }

    // test 2D problem
    /// \todo - This test fails with current tolerance.
    void xTestSimpleDg0ParabolicAssembler2DZeroDirichWithSourceTermOnFineMeshWithSmallDt()
    {
        // Create mesh from mesh reader
        FemlabMeshReader<2,2> mesh_reader("mesh/test/data/",
                                          "femlab_square_nodes.dat",
                                          "femlab_square_elements.dat",
                                          "femlab_square_edges.dat");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquationWithSourceTerm<2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorEnd();

        while (iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            ConstBoundaryCondition<2> *p_dirichlet_boundary_condition =
                new ConstBoundaryCondition<2>(-0.25*(x*x+y*y));
            bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            iter++;
        }

        // Assembler
        SimpleDg0ParabolicAssembler<2,2,true> assembler(&mesh,&pde,&bcc);

        // initial condition, u(0,x) = sin(x*pi)*sin(y*pi)-0.25*(x^2+y^2);
        Vec initial_condition = PetscTools::CreateVec(mesh.GetNumNodes());

        double *p_initial_condition;
        VecGetArray(initial_condition, &p_initial_condition);

        int lo, hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNode(global_index)->GetPoint()[0];
            double y = mesh.GetNode(global_index)->GetPoint()[1];
            p_initial_condition[local_index] = sin(x*M_PI)*sin(y*M_PI)-0.25*(x*x+y*y);
        }
        VecRestoreArray(initial_condition, &p_initial_condition);

        double t_end = 0.1;
        assembler.SetTimes(0, t_end, 0.001);
        assembler.SetInitialCondition(initial_condition);

        Vec result = assembler.Solve();

        // Check result
        double *p_result;
        VecGetArray(result, &p_result);

        // Solution should be u = e^{-t*2*pi*pi} sin(x*pi) sin(y*pi) - 0.25(x^2+y^2), t=0.1
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNode(global_index)->GetPoint()[0];
            double y = mesh.GetNode(global_index)->GetPoint()[1];
            double u = exp(-0.1*2*M_PI*M_PI)*sin(x*M_PI)*sin(y*M_PI)-0.25*(x*x+y*y);
            TS_ASSERT_DELTA(p_result[local_index], u, 0.001);
        }
        VecRestoreArray(result, &p_result);
        VecDestroy(initial_condition);
        VecDestroy(result);
    }

    // test 2D problem
    void TestSimpleDg0ParabolicAssembler2DNeumannOnCoarseMesh()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();

        while (iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];

            if ((fabs(y) < 0.01) || (fabs(y - 1.0) < 0.01) || (fabs(x) < 0.01))
            {
                ConstBoundaryCondition<2> *p_dirichlet_boundary_condition
                = new ConstBoundaryCondition<2>(x);
                bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            }

            iter++;
        }

        TetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2> *p_neumann_boundary_condition =
            new ConstBoundaryCondition<2>(1.0);

        while (surf_iter < mesh.GetBoundaryElementIteratorEnd())
        {
            int node = (*surf_iter)->GetNodeGlobalIndex(0);
            double x = mesh.GetNode(node)->GetPoint()[0];

            if (fabs(x - 1.0) < 0.01)
            {
                bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
            }

            surf_iter++;
        }

        // Assembler
        SimpleDg0ParabolicAssembler<2,2,true> assembler(&mesh,&pde,&bcc);

        // initial condition, u(0,x,y) = sin(0.5*M_PI*x)*sin(M_PI*y)+x
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(0.5*M_PI*x)*sin(M_PI*y)+x;;
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);


        double t_end = 0.1;
        assembler.SetTimes(0, t_end, 0.01);
        assembler.SetInitialCondition(initial_condition);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // check result
        // Solution should be u = e^{-5/4*M_PI*M_PI*t} sin(0.5*M_PI*x)*sin(M_PI*y)+x, t=0.1
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = exp((-5/4)*M_PI*M_PI*0.1) * sin(0.5*M_PI*x) * sin(M_PI*y) + x;
            TS_ASSERT_DELTA(result_repl[i], u, u*0.15);
        }

        VecDestroy(initial_condition);
        VecDestroy(result);
    }

    // test 2D problem
    void TestSimpleDg0ParabolicAssembler2DNeumann()
    {
        // Create mesh from mesh reader
        FemlabMeshReader<2,2> mesh_reader("mesh/test/data/",
                                          "femlab_square_nodes.dat",
                                          "femlab_square_elements.dat",
                                          "femlab_square_edges.dat");

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();

        while (iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];

            if ((fabs(y) < 0.01) || (fabs(y - 1.0) < 0.01) || (fabs(x) < 0.01))
            {
                ConstBoundaryCondition<2> *p_dirichlet_boundary_condition =
                    new ConstBoundaryCondition<2>(x);
                bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            }

            iter++;
        }

        TetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2> *p_neumann_boundary_condition =
            new ConstBoundaryCondition<2>(1.0);

        while (surf_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            int node = (*surf_iter)->GetNodeGlobalIndex(0);
            double x = mesh.GetNode(node)->GetPoint()[0];

            if (fabs(x - 1.0) < 0.01)
            {
                bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
            }

            surf_iter++;
        }

        // Assembler
        SimpleDg0ParabolicAssembler<2,2,true> assembler(&mesh,&pde,&bcc);

        // initial condition, u(0,x,y) = sin(0.5*M_PI*x)*sin(M_PI*y)+x
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(0.5*M_PI*x)*sin(M_PI*y)+x;
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);

        assembler.SetTimes(0, 0.1, 0.01);
        assembler.SetInitialCondition(initial_condition);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // check result
        // Solution should be u = e^{-5/4*M_PI*M_PI*t} sin(0.5*M_PI*x)*sin(M_PI*y)+x, t=0.1
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double u = exp((-5/4)*M_PI*M_PI*0.1) * sin(0.5*M_PI*x) * sin(M_PI*y) + x;
            TS_ASSERT_DELTA(result_repl[i], u, u*0.1);
        }

        VecDestroy(initial_condition);
        VecDestroy(result);
    }

    void TestHeatEquationSolutionDoesntDrift2D()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        // Boundary conditions - non-zero constant dirichlet on boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        ConstBoundaryCondition<2>* dirichlet_bc = new ConstBoundaryCondition<2>(-84.5);
        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            bcc.AddDirichletBoundaryCondition(*iter, dirichlet_bc);
            iter++;
        }

        // Assembler
        SimpleDg0ParabolicAssembler<2,2,true> assembler(&mesh,&pde,&bcc);

        // initial condition;
        Vec initial_condition = PetscTools::CreateVec(mesh.GetNumNodes(), -84.5);

        double t_end = 1.0;
        assembler.SetTimes(0, t_end, 0.01);
        assembler.SetInitialCondition(initial_condition);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // Check solution is constant throughout the mesh
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            TS_ASSERT_DELTA(result_repl[i],-84.5, 0.0002);
        }

        VecDestroy(initial_condition);
        VecDestroy(result);
    }

    void TestHeatEquationSolutionDoesntDrift1D()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<1> pde;

        // Boundary conditions - non-zero constant dirichlet on boundary
        BoundaryConditionsContainer<1,1,1> bcc;
        TetrahedralMesh<1,1>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        ConstBoundaryCondition<1>* dirichlet_bc = new ConstBoundaryCondition<1>(-84.5);
        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            bcc.AddDirichletBoundaryCondition(*iter, dirichlet_bc);
            iter++;
        }

        // Assembler - created using Set methods for coverage
        SimpleDg0ParabolicAssembler<1,1,true> assembler(NULL,&pde,NULL);
        assembler.SetMesh(&mesh);
        assembler.SetBoundaryConditionsContainer(&bcc);

        // initial condition;
        Vec initial_condition = PetscTools::CreateVec(mesh.GetNumNodes(), -84.5);

        double t_end = 1;
        assembler.SetTimes(0, t_end, 0.01);
        assembler.SetInitialCondition(initial_condition);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // Check solution is constant throughout the mesh
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            TS_ASSERT_DELTA(result_repl[i],-84.5, 0.0001);
        }

        VecDestroy(initial_condition);
        VecDestroy(result);
    }

    // commented out heat equation with 2d mesh and initial condition non-zero at centre,
    // writing out data (doesn't test anything, wanted to see if we get a circular
    // diffusion pattern on such a small mesh, to compare with monodomain with
    // centre stimulus - result doesn't look like a circle)
    // !Need to change the diffusion coefficient to 0.001 if running this!
    void DONOT_TestSimpleDg0ParabolicAssembler2DZeroNeumannNonZeroInCentre()
    {
        // read mesh on [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        HeatEquation<2> pde;

        BoundaryConditionsContainer<2,2,1> bcc;
        TetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2> *p_neumann_boundary_condition =
            new ConstBoundaryCondition<2>(0.0);

        while (surf_iter < mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
            surf_iter++;
        }

        // Assembler
        SimpleDg0ParabolicAssembler<2,2,true> assembler(&mesh,&pde,&bcc);

        // initial condition;
        // choose initial condition sin(x*pi)*sin(y*pi) as this is an eigenfunction of
        // the heat equation.
        Vec initial_condition = PetscTools::CreateVec(mesh.GetNumNodes());

        double *p_initial_condition;
        VecGetArray(initial_condition, &p_initial_condition);

        int lo, hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);

        // stimulate
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            p_initial_condition[local_index] = 0;
            switch (global_index)
            {
                case 60:
                case 165:
                case 166:
                case 175:
                case 176:
                    p_initial_condition[local_index] = 100;
                    break;
            }
        }
        VecRestoreArray(initial_condition, &p_initial_condition);

        double time = 0;
        double t_end = 0.1;
        double dt = 0.001;
        assembler.SetInitialCondition(initial_condition);

        int time_var_id = 0;
        int heat_var_id = 0;

        ParallelColumnDataWriter *p_test_writer;
        p_test_writer = new ParallelColumnDataWriter("2DHeatEquation", "2DHeatEquation");

        p_test_writer->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
        time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");

        heat_var_id = p_test_writer->DefineVariable("T","K");
        p_test_writer->EndDefineMode();

        p_test_writer->PutVariable(time_var_id, time);
        p_test_writer->PutVector(heat_var_id, initial_condition);
        p_test_writer->AdvanceAlongUnlimitedDimension();

        Vec result;

        while (time < t_end)
        {
            time += dt;
            assembler.SetTimes(time, time+dt, dt);

            result = assembler.Solve();

            assembler.SetInitialCondition(result);

            p_test_writer->PutVariable(time_var_id, time);
            p_test_writer->PutVector(heat_var_id, result);
            p_test_writer->AdvanceAlongUnlimitedDimension();
        }
        VecDestroy(initial_condition);
        VecDestroy(result);
    }
};

#endif //_TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_
