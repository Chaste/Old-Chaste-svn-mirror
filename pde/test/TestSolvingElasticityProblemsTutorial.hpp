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
#ifndef TESTSOLVINGELASTICITYPROBLEMS_HPP_
#define TESTSOLVINGELASTICITYPROBLEMS_HPP_

/*
 * = Solving solid mechanics problems =
 *
 * In this tutorial we show how Chaste can be used to solve solid mechanics problems.
 * We assume some the reader has some familiarity with solid mechanics problems
 * (note: the equations of nonlinear elasticity are given in Section 4.1 of the PDF on equations
 * and finite element implementations in ChasteGuides -> Miscellaneous information). In brief,
 * there several facets to such models:
 *  * Time-dependent problems versus static problems
 *  * Linear elasticity versus nonlinear elasticity
 *  * Compressible versus incompressible materials
 *  * The type of material behaviour (elastic, visco-elastic, etc..)
 *  * Specification of geometry, material law, body force, displacement boundary conditions, traction boundary conditions
 *
 * EMPTYLINE
 *
 * The solvers currently implemented are STATIC (time-independent) and use NONLINEAR ELASTICITY. The main solver
 * solves for an INCOMPRESSIBLE deformation, although there is now a COMPRESSIBLE solver. The material behaviour is
 * assumed to be ELASTIC (stress is just a function of strain, not strain-rate etc), and in particular HYPER-ELASTIC
 * (stress is a function of strain via a 'strain energy function' (S.E.F.), for which stress is obtained by differentiating the
 * S.E.F. with respect to strain).
 *
 * EMPTYLINE
 *
 * To solve a mechanics problem we need to
 *  * Choose the solver (compressible or incompressible)
 *  * Specify the geometry (ie the mesh)
 *  * Specify the material law (ie the strain-energy function)
 *  * Specify the BODY FORCE -- this is a force (per unit volume) acting throughout the body (ie, gravity)
 *  * Specify some DISPLACEMENT BOUNDARY CONDITIONS -- some part of the boundary must have the displacement specified on it
 *  * Specify TRACTION BOUNDARY CONDITIONS (if non-zero) on the rest of the boundary -- tractions are pressures applied
 *  the rest of the surface of the deformable object.
 *
 * EMPTYLINE
 *
 *  '''VERY IMPORTANT NOTE:''' Make sure you read the comment about HYPRE below before going to 3D or refining the meshes in these tests.
 *
 * EMPTYLINE
 *
 * As always we include this first class as a test suite */
#include <cxxtest/TestSuite.h>
/* On some systems there is a clash between Boost Ublas includes and PETSc.  This can be
 * resolved by making sure that Chaste's interface to the Boost libraries are included
 * as early as possible.
 */
#include "UblasCustomFunctions.hpp"
/* The incompressible solver is called `NonlinearElasticitySolver` */
#include "NonlinearElasticitySolver.hpp"
/* The compressible solver is called `CompressibleNonlinearElasticitySolver` */
#include "CompressibleNonlinearElasticitySolver.hpp"
/* The simplest incompressible material law is the Mooney-Rivlin material law (of which
 * Neo-Hookean laws are a subset)  */
#include "MooneyRivlinMaterialLaw.hpp"
/* Another incompressible material law */
#include "ExponentialMaterialLaw.hpp"
/* This is a useful helper class */
#include "NonlinearElasticityTools.hpp"
/* As before: !PetscSetupAndFinalize.hpp must be included in every test that uses PETSc. Note that it
 * cannot be included in the source code. */
#include "PetscSetupAndFinalize.hpp"


/* Ignore this function until later in the tutorial */
c_vector<double,2> MyTraction(c_vector<double,2>& X, double time)
{
    c_vector<double,2> traction = zero_vector<double>(2);
    traction(0) = X(0);
    return traction;
}



/*
 *
 * EMPTYLINE
 *
 * == Simple incompressible deformation: 2D shape hanging under gravity ==
 *
 * EMPTYLINE
 */
class TestSolvingElasticityProblemsTutorial : public CxxTest::TestSuite
{
public:
    /* In the first test we use INCOMPRESSIBLE nonlinear elasticity. For such problems there is a constraint
     * on the deformation, which results in a pressure field (a Lagrange multiplier) which needs to be solved
     * for together with the deformation.
     *
     * EMPTYLINE
     *
     * All the mechanics solvers solve for the deformation using the finite element method with QUADRATIC
     * basis functions for the deformation: this necessitates the use of a QUADRATIC MESH (such meshes have
     * extra nodes that aren't vertices of elements, in this case midway along each edge). The displacement
     * is solved for at ''each node'' in the mesh (including internal [non-vertex] nodes), whereas the pressure
     * is only solved for at each vertex. (In FEM terms, quadratic interpolation for displacement, linear
     * interpolation for pressure, which is required for stability).
     *
     * Note: 1D incompressible solves are meaningless and therefore not allowed.
     *
     */
    void TestSimpleIncompressibleProblem() throw(Exception)
    {
        /* First, define the geometry. This should be specified using the `QuadraticMesh` class, which inherits from `TetrahedralMesh`
         * and has mostly the same interface. Here we define a 0.8 by 1 rectangle, with elements 0.1 wide.
         * (`QuadraticMesh`s can also be read in using `TrianglesMeshReader`; see rest of code base for examples of this).
         */
        QuadraticMesh<2> mesh;
        mesh.ConstructRegularSlabMesh(0.1 /*stepsize*/, 0.8 /*width*/, 1.0 /*height*/);

        /* We use a Mooney-Rivlin material law, which applies to isotropic materials and has two parameters.
         * Restricted to 2d however, it only has one parameter, which can be thought of as the total
         * stiffness. We declare a Mooney-Rivlin law, setting the parameter to 1.
         */
        MooneyRivlinMaterialLaw<2> law(1.0);

        /* Next, the body force (per unit volume). In realistic problems this will either be
         * acceleration due to gravity (ie b=(0,-9.81)) or zero if the effect of gravity can be neglected.
         * In this problem we apply a gravity-like downward force.
         */
        c_vector<double,2> body_force;
        body_force(0) =  0.0;
        body_force(1) = -2.0;

        /* The density of the body is also needed */
        double density = 1.0;

        /* Two types of boundary condition are required: displacement and traction. As with the other PDE solvers,
         * the displacement (Dirichlet) boundary conditions are specified at nodes, whereas traction (Neumann) boundary
         * conditions are specified on boundary elements.
         *
         * EMPTYLINE
         *
         * In this test we apply displacement boundary conditions on one surface of the mesh, the upper (Y=1.0) surface.
         * We do not specify what the displacement is, which means zero-displacement will be prescribed for these nodes.
         * We do not specify any traction boundary conditions, which means that (effectively) zero-traction boundary
         * conditions (ie zero pressures) are applied on the three other surfaces.
         *
         * EMPTYLINE
         *
         * We need to get a `std::vector` of all the node indices that we want to fix. The `NonlinearElasticityTools`
         * has a static method for helping do this: the following gets all the nodes for which Y=1.0. The second
         * argument (the '1') indicates Y (eg, `GetNodesByComponentValue(mesh, 0, 10)` would correspond to X=10).
         */
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh, 1, 1.0);

        /* Now we create the (incompressible) solver, passing in the mesh, law, force, density, output directory
         * name and fixed nodes..
         */
        NonlinearElasticitySolver<2> solver(&mesh,
                                            &law,
                                            body_force,
                                            density,
                                            "SimpleIncompressibleElasticityTutorial",
                                            fixed_nodes);

        /* .. and call `Solve()` */
        solver.Solve();

        /* This test is just here to (help) check nothing has gotten changed in this test. Note that since we are solving
         * a nonlinear problem we have to a nonlinear solver. We use Newton's method (with damping). In this test
         * 4 iterations were needed to converge. */
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 4u);

        /* '''Visualisation'''. Go to the folder `SimpleIncompressibleElasticityTutorial` in your test-output directory.
         * There should be 5 files, named solution_[i].nodes, for i=0,1,2,3,4. These are the solutions after each Newton
         * iteration. solution_0.nodes is therefore the initial condition, and solution_4 is the final solution. Each file
         * has two columns, the x and y locations of each node. To visualise the solution in say Matlab or Octave, you could do:
         * `x=load('solution_4.nodes'); plot(x(:,1),x(:,2),'k*')`.
         *
         * EMPTYLINE
         *
         * To get the actual solution from the solver, use these two methods. Note that the first
         * gets the deformed position (ie the new location, not the displacement), and will be of size
         * num_total_nodes; the second will be of size num_vertices.
         */
        std::vector<c_vector<double,2> >& r_deformed_positions = solver.rGetDeformedPosition();
        std::vector<double>& r_pressures = solver.rGetPressures();
        /* Let us obtain the values of the new position, and the pressure, at the bottom right corner node. */
        unsigned node_index = 8;
        assert( fabs(mesh.GetNode(node_index)->rGetLocation()[0] - 0.8) < 1e-6); // check that X=0.8, ie that we have the correct node,
        assert( fabs(mesh.GetNode(node_index)->rGetLocation()[1] - 0.0) < 1e-6); // check that Y=0.0, ie that we have the correct node,
        std::cout << "New position: " << r_deformed_positions[node_index](0) << " " << r_deformed_positions[node_index](1) << "\n";
        std::cout << "Pressure: " << r_pressures[node_index] << "\n";



//
////TODO - implement this inside solver and add call here
//        CmguiDeformedSolutionsWriter<2> writer("SimpleIncompressibleElasticityTutorial/cmgui",
//                                                "solution",
//                                                mesh,
//                                                WRITE_QUADRATIC_MESH);
//
//        // writes solution_0.exnode and solution_0.exelem using the mesh
//        writer.WriteInitialMesh();
//        writer.WriteDeformationPositions(r_deformed_positions, 1);
//        writer.WriteDeformationPositions(r_deformed_positions, 2);
//        writer.WriteDeformationPositions(r_deformed_positions, 3);
//        writer.WriteDeformationPositions(r_deformed_positions, 4);
//        writer.WriteCmguiScript();
    }

    /*
     * EMPTYLINE
     *
     * == Incompressible deformation: 2D shape hanging under gravity with balancing traction ==
     *
     * EMPTYLINE
     *
     * We now repeat the above test but include a traction, on the bottom surface (Y=0). We apply this
     * in the inward direction so that is counters (somewhat) the effect of gravity.
     */
    void TestIncompressibleProblemWithTractions() throw(Exception)
    {
        /* All of this is exactly as above */
        QuadraticMesh<2> mesh;
        mesh.ConstructRegularSlabMesh(0.1 /*stepsize*/, 0.8 /*width*/, 1.0 /*height*/);

        MooneyRivlinMaterialLaw<2> law(1.0);

        c_vector<double,2> body_force;
        body_force(0) =  0.0;
        body_force(1) = -2.0;

        double density = 1.0;

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh, 1, 1.0);

        /* Now the traction boundary conditions. We need to collect all the boundary elements on the surface which we want to
         * apply non-zero tractions, put them in a `std::vector`, and create a corresponding `std::vector` of the tractions
         * for each of the boundary elements. Note that the each traction is a 2D vector with dimensions of pressure.
         *
         * EMPTYLINE
         *
         * Declare the data structures
         */
        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        /* Create a constant traction */
        c_vector<double,2> traction;
        traction(0) = 0;
        traction(1) = 1.0; // this choice of sign corresponds to an inward force (if applied to the bottom surface)
        /* Loop over boundary elements */
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
                = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            /* If the centre of the element has Y value of 0.0, it is on the surface we need */
            if (fabs((*iter)->CalculateCentroid()[1] - 0.0)<1e-6)
            {
                /* Put the boundary element and the constant traction into the stores. */
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }
        /* A quick check */
        assert(boundary_elems.size()==8u);
        /* Create solver as before */
        NonlinearElasticitySolver<2> solver(&mesh,
                                            &law,
                                            body_force,
                                            density,
                                            "IncompressibleElasticityWithTractionsTutorial",
                                            fixed_nodes);
        /* This is how to set the tractions */
        solver.SetSurfaceTractionBoundaryConditions(boundary_elems, tractions);
        /* Call `Solve()` */
        solver.Solve();

        /* Another quick check */
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 3u); // 3 rather than 4 this time

        /* Visualise as before by going to `IncompressibleElasticityWithTractionsTutorial` and doing
         * `x=load('solution_3.nodes'); plot(x(:,1),x(:,2),'m*')`. The effect of the traction should be
         * clear (especially when compared to the results of the first test).
         */
    }

    /* == Incompressible deformation: non-zero displacement boundary conditions, functional tractions ==
     *
     * We now consider a more complicated example. We prescribe particular new locations for the nodes
     * on the Dirichlet boundary, and also show how to prescribe a traction that is given in functional form
     * rather than prescribed for each boundary element.
     */
    void TestIncompressibleProblemMoreComplicatedExample() throw(Exception)
    {
        /* Create a mesh */
        QuadraticMesh<2> mesh;
        mesh.ConstructRegularSlabMesh(0.1 /*stepsize*/, 1.0 /*width*/, 1.0 /*height*/);

        /* Use a different material law this time, an exponential material law.
         * The material law needs to inherit from `AbstractIncompressibleMaterialLaw`,
         * and there are a few implemented, see `pde/src/problem` */
        ExponentialMaterialLaw<2> law(1.0, 0.5); // First parameter is 'a', second 'b', in W=a*exp(b(I1-3))
        /* Note that it is possible to specify different material laws for each element in the mesh (for example
         * for using different stiffnesses in different regions. Just create a `std::vector<AbstractIncompressibleMaterial<DIM>*>`
         * and fill it in with the material law for each element, and pass as the first argument of the solver.
         *
         * EMPTYLINE
         *
         * Now specify the fixed nodes, and their new locations. Create `std::vector`s for each. */
        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > locations;
        /* Loop over the mesh nodes */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {

            /* If the node is on the Y=0 surface (the LHS) */
            if ( fabs(mesh.GetNode(i)->rGetLocation()[1])<1e-6)
            {
                /* Add it to the list of fixed nodes */
                fixed_nodes.push_back(i);
                /* and define a new position x=(X,0.1*X^2^) */
                c_vector<double,2> new_location;
                double X = mesh.GetNode(i)->rGetLocation()[0];
                new_location(0) = X;
                new_location(1) = 0.1*X*X;
                locations.push_back(new_location);
            }
        }

        /* Now collect all the boundary elements on the top surface, as before, except
         * here we don't create the tractions for each element
         */
        std::vector<BoundaryElement<1,2>*> boundary_elems;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
                = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            /* If Y=1, have found a boundary element */
            if (fabs((*iter)->CalculateCentroid()[1] - 1.0)<1e-6)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
            }
        }

        /* Create the solver, passing in the new locations as the final argument */
        NonlinearElasticitySolver<2> solver(&mesh,
                                            &law,
                                            zero_vector<double>(2),
                                            1.0,
                                            "IncompressibleElasticityMoreComplicatedExample",
                                            fixed_nodes,
                                            &locations);

        /* Now call `SetFunctionalTractionBoundaryCondition`, passing in the boundary elements,
         * together with a ''function pointer'' (just the name of the function) to a
         * function returning traction in terms of position (and time [see below]). This function is defined
         * above, before the tests. It has take in a `c_vector` (position) and a double (time), and return a
         * `c_vector` (traction), and will only be called using the quadrature points in the boundary elements being
         * passed in. The function `MyTraction` above defines a horizontal traction (ie a shear stress, since it is
         * applied to the top surface) which increases in magnitude across the object.
         */
        solver.SetFunctionalTractionBoundaryCondition(boundary_elems, MyTraction);
        /* Note: You can also call `solver.SetFunctionalBodyForce(MyBodyForce)`, though is only
         * really useful for constructing problems with exact solutions.
         *
         * EMPTYLINE
         *
         * Call `Solve()` */
        solver.Solve();

        /* Another quick check */
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 6u);
        /* Visualise with `x=load('solution_6.nodes'); plot(x(:,1),x(:,2),'b*')`
         *
         * EMPTYLINE
         *
         * '''Advanced:''' Note that the function `MyTraction` takes in time, which it didn't use. In the above it would have been called
         * with t=0. The current time can be set using `SetCurrentTime()`. The idea is that the user may want to solve a
         * sequence of static problems with time-dependent tractions (say), for which they should allow `MyTraction` to
         * depend on time, and put the solve inside a time-loop, for example:
         */
         //for(double t=0; t<T; t+=dt)
         //{
         //    solver.SetCurrentTime(t);
         //    solver.Solve();
         //}
         /* In this the current time would be passed through to `MyTraction`.*/
    }

    /* == IMPORTANT: Using HYPRE ==
     *
     * When running problems in 3D, or with more elements, it is vital to change the linear solver to use HYPRE, an algebraic multigrid
     * solver. Without HYRPE, the linear solve (i) may become very very slow; or (ii) may not converge, in which case the nonlinear
     * solve will (probably) not converge. HYPRE is currently not a pre-requisite for installing Chaste, hence this is not (currently)
     * the default linear solver for mechanics problems. You need to have Petsc installed with HYPRE. However, if you followed installation
     * instructions for Chaste 2.1 or later, you probably do already have Petsc installed with HYPRE.
     *
     * EMPTYLINE
     *
     * To switch on HYPRE, open the file `pde/src/solver/AbstractNonlinearElasticitySolver` and uncomment the line
     * #define MECH_USE_HYPRE
     * near the top of the file (currently: line 53).
     *
     * EMPTYLINE
     *
     * Mechanics solves being nonlinear are ''expensive'', so it is recommended you also use `build=GccOpt_ndebug` (when running scons)
     * on larger problems.
     *
     * EMPTYLINE
     *
     * Note: Petsc unfortunately doesn't quit if you try to use HYPRE without it being installed, but it spew lots of error messages.
     *
     * EMPTYLINE
     *
     * == Compressible Problems ==
     *
     * EMPTYLINE
     *
     * To solve compressible elasticity problems, all that needs to be changed is to use `CompressibleNonlinearElasticitySolver` instead
     * of `NonlinearElasticitySolver` (making sure we include it), and noting that there is no pressure computed. See
     * `TestCompressibleNonlinearElasticitySolver`. Compressible solid mechanics is in the process of being implemented properly, tutorials
     * to be added later.
     */
};



#endif /*TESTSOLVINGELASTICITYPROBLEMS_HPP_*/
