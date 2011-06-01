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

#ifndef TESTLINEARPARABOLICPDESYSTEMWITHCOUPLEDODESYSTEMSOLVER_HPP_
#define TESTLINEARPARABOLICPDESYSTEMWITHCOUPLEDODESYSTEMSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainerImplementation.hpp"
#include "AbstractBoundaryConditionsContainerImplementation.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OdeSystemForCoupledHeatEquation.hpp"
#include "OdeSystemForCoupledHeatEquationWithSource.hpp"
#include "HeatEquationForCoupledOdeSystem.hpp"
#include "HeatEquationWithSourceForCoupledOdeSystem.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"

class TestLinearParabolicPdeSystemWithCoupledOdeSystemSolver : public CxxTest::TestSuite
{
public:

    void TestHeatEquationWithSourceWithCoupledOdeSystemIn1dWithZeroDirichlet()
    {
        // Create mesh of domain [0,1]
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create PDE system object
        HeatEquationWithSourceForCoupledOdeSystem<1> pde;

        // Define zero Neumann boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
        iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        // Create the correct number of ODE systems
        double a = 5.0;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> ode_systems;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ode_systems.push_back(new OdeSystemForCoupledHeatEquationWithSource(a));
        }

        // Create PDE system solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<1,1,1> solver(&mesh, &pde, &bcc, ode_systems);

        // Test setting end time and timestep
        TS_ASSERT_THROWS_THIS(solver.SetTimes(1.0, 0.0), "Starting time has to less than ending time");
        TS_ASSERT_THROWS_THIS(solver.SetTimeStep(0.0), "Time step has to be greater than zero");

        // Set end time and timestep
        double t_end = 0.5;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.001);

        // Set initial condition u(x,0) = 1 + cos(pi*x)
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            init_cond[i] = 1 + cos(M_PI*x);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE system and store result
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        /*
         * Test that solution is given by
         *
         * u(x,t) = 1 + (1 - exp(-a*t))/a + exp(-pi*pi*t)*cos(pi*x),
         * v(x,t) = exp(-a*t),
         *
         * with t = t_end.
         */
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = 1 + (1 - exp(-a*t_end))/a + exp(-M_PI*M_PI*t_end)*cos(M_PI*x);
            TS_ASSERT_DELTA(result_repl[i], u, 0.1);

            double v = exp(-a*t_end);
            TS_ASSERT_DELTA(ode_systems[i]->rGetStateVariables()[0], v, 0.1);
        }

        // Tidy up
        VecDestroy(initial_condition);
        VecDestroy(result);
    }

    void TestHeatEquationWithCoupledOdeSystemIn2dWithZeroDirichlet()
    {
        // Create mesh of the domain [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create PDE system object
        HeatEquationForCoupledOdeSystem<2> pde;

        // Define zero Dirichlet boundary conditions on entire boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Create the correct number of ODE systems
        double a = 5.0;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> ode_systems;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ode_systems.push_back(new OdeSystemForCoupledHeatEquation(a));
        }

        // Create PDE system solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,1> solver(&mesh, &pde, &bcc, ode_systems);

        // Set end time and timestep
        double t_end = 0.5;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.001);

        /*
         * Set initial condition
         *
         * u(x,y,0) = sin(pi*x)*sin(pi*y),
         *
         * which is an eigenfunction of the heat equation.
         */
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(M_PI*x)*sin(M_PI*y);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE system and store result
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        /*
         * Test that solution is given by
         *
         * u(x,y,t) = e^{-2*pi*pi*t} sin(pi*x)*sin(pi*y),
         * v(x,y,t) = 1 + (1 - e^{-2*pi*pi*t})*sin(pi*x)*sin(pi*y)*a/(2*pi*pi),
         *
         * with t = t_end.
         */
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];

            double u = exp(-2*M_PI*M_PI*t_end)*sin(M_PI*x)*sin(M_PI*y);
            double v = 1.0 + (a/(2*M_PI*M_PI))*(1 - exp(-2*M_PI*M_PI*t_end))*sin(M_PI*x)*sin(M_PI*y);

            TS_ASSERT_DELTA(result_repl[i], u, 0.01);
            TS_ASSERT_DELTA(ode_systems[i]->rGetStateVariables()[0], v, 0.01);
        }

        // Tidy up
        VecDestroy(initial_condition);
        VecDestroy(result);
    }

    void TestSolveAndWriteResultsToFileMethod()
    {
        // Create mesh of the domain [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create PDE system object
        HeatEquationForCoupledOdeSystem<2> pde;

        // Define zero Dirichlet boundary conditions on entire boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Create the correct number of ODE systems
        double a = 5.0;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> ode_systems;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ode_systems.push_back(new OdeSystemForCoupledHeatEquation(a));
        }

        // Create PDE system solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,1> solver(&mesh, &pde, &bcc, ode_systems);

        // Set output directory
        solver.SetOutputDirectory("TestHeatEquationForCoupledOdeSystemIn2dWithZeroDirichletWithOutput");

        // Set end time and timestep
        double t_end = 1.0;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.01);

        // Set sampling timestep
        solver.SetSamplingTimeStep(0.1);

        /*
         * Set initial condition
         *
         * u(x,y,0) = sin(pi*x)*sin(pi*y),
         *
         * which is an eigenfunction of the heat equation.
         */
        std::vector<double> init_cond(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            init_cond[i] = sin(M_PI*x)*sin(M_PI*y);
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetInitialCondition(initial_condition);

        TS_ASSERT_THROWS_NOTHING(solver.SolveAndWriteResultsToFile());

        // Tidy up
        VecDestroy(initial_condition);
    }
};

#endif /*TESTLINEARPARABOLICPDESYSTEMWITHCOUPLEDODESYSTEMSOLVER_HPP_*/
