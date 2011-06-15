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
/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTSOLVINGLINEARPARABOLICPDESYSTEMSWITHCOUPLEDODESYSTEMSTUTORIAL_HPP_
#define TESTSOLVINGLINEARPARABOLICPDESYSTEMSWITHCOUPLEDODESYSTEMSTUTORIAL_HPP_

/*
 * = Examples showing how to solve a system of coupled linear parabolic PDEs and ODEs =
 *
 * In this tutorial we show how Chaste can be used to solve a system of coupled linear
 * parabolic PDEs and ODEs. This test uses the {{{LinearParabolicPdeSystemWithCoupledOdeSystemSolver}}}.
 *
 * EMPTYLINE
 *
 * The following header files need to be included.
 * First we include the header needed to define this class as a test suite.
 */
#include <cxxtest/TestSuite.h>
/*
 * On some systems there is a clash between Boost Ublas includes and PETSc.  This can be
 * resolved by making sure that Chaste's interface to the Boost libraries are included
 * as early as possible.
 */
#include "UblasCustomFunctions.hpp"
/*
 * This is the class that is needed to solve a system of coupled linear
 * parabolic PDEs and ODEs.
 */
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
/*
 * The next header file defines the Schnackenberg system, which comprises
 * two reaction-diffusion PDEs that are coupled through their reaction terms.
 */
#include "SchnackenbergCoupledPdeSystem.hpp"
/*
 * The next header file will allow us to specify a random initial condition.
 */
#include "RandomNumberGenerator.hpp"
/*
 * We then include header files that allow us to specify boundary conditions for the PDEs,
 * deal with meshes and output files, and use PETSc. As noted before, !PetscSetupAndFinalize.hpp
 * must be included in every test that uses PETSc.
 */
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"

/* == Test 1: Solving the Schnackenberg system ==
 *
 * Here, we solve the Schnackenberg system...
 *
 * EMPTYLINE
 *
 * Next, we define the test suite (a class). It is sensible to name it the same
 * as the filename. The class should inherit from {{{CxxTest::TestSuite}}}.
 */
class TestSolvingLinearParabolicPdeSystemsWithCoupledOdeSystemsTutorial : public CxxTest::TestSuite
{
/*
 * All individual test defined in this test suite '''must''' be declared as public.
 */
public:
    /*
     * Define a particular test.
     */
    void TestSchnackenbergCoupledPdeSystemIn1dWithNonZeroDirichlet()
    {
        /*
         * First we declare a mesh reader which reads mesh data files of the 'Triangle'
         * format. The path given is the relative to the main Chaste directory. The reader
         * will look for three datafiles, [name].nodes, [name].ele and (in 2d or 3d)
         * [name].edge. Note that the first template argument here is the dimension of the
         * elements in the mesh ({{{ELEMENT_DIM}}}), and the second is the dimension of the nodes,
         * i.e. the dimension of the space the mesh lives in ({{{SPACE_DIM}}}). Usually
         * {{{ELEMENT_DIM}}} and {{{SPACE_DIM}}} will be equal. */
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_1000_elements");
        /* Now declare a tetrahedral mesh with the same dimensions */
        TetrahedralMesh<1,1> mesh;
        /* Construct the mesh using the mesh reader */
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create PDE system object
        SchnackenbergCoupledPdeSystem<1> pde(1e-4, 1e-2, 0.1, 0.2, 0.3, 0.1);

        /* A set of boundary conditions are stored in a {{{BoundaryConditionsContainer}}}. The
         * three template arguments are ELEMENT_DIM, SPACE_DIM and PROBLEM_DIM, the latter being
         * the number of unknowns we are solving for. We have two unknown (ie u is a 2-vector),
         * so in this case {{{PROBLEM_DIM}}}=2. */
        BoundaryConditionsContainer<1,1,2> bcc;
        ConstBoundaryCondition<1>* p_bc_for_u = new ConstBoundaryCondition<1>(2.0);
        ConstBoundaryCondition<1>* p_bc_for_v = new ConstBoundaryCondition<1>(0.75);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_bc_for_u, 0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_bc_for_v, 1);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(mesh.GetNumNodes()-1), p_bc_for_u, 0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(mesh.GetNumNodes()-1), p_bc_for_v, 1);

        // Create PDE system solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<1,1,2> solver(&mesh, &pde, &bcc);

        // Set end time and time step
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-1);

        // Create initial conditions that are random perturbations of the uniform steady state
        std::vector<double> init_conds(2*mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            init_conds[2*i] = fabs(2.0 + RandomNumberGenerator::Instance()->ranf());
            init_conds[2*i + 1] = fabs(0.75 + RandomNumberGenerator::Instance()->ranf());
        }
        Vec initial_condition = PetscTools::CreateVec(init_conds);
        solver.SetInitialCondition(initial_condition);

        // Solve PDE system and store result
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        // Store results in an accessible form
        ReplicatableVector solution_repl(result);

        /*
         * All Petsc {{{Vec}}}s should be destroyed when they are no longer needed.
         */
        VecDestroy(initial_condition);
        VecDestroy(result);
    }
};
#endif /*TESTSOLVINGLINEARPARABOLICPDESYSTEMSWITHCOUPLEDODESYSTEMSTUTORIAL_HPP_*/
