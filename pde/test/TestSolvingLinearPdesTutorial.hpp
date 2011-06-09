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
#ifndef TESTSOLVINGLINEARPDESTUTORIAL_HPP_
#define TESTSOLVINGLINEARPDESTUTORIAL_HPP_

/*
 * = Examples showing how to solve linear elliptic and parabolic PDEs =
 *
 * In this tutorial we show how Chaste can be used to solve linear PDEs. The first test
 * uses the {{{SimpleLinearEllipticSolver}}} to solve a linear elliptic PDE, and the
 * second test uses the {{{SimpleLinearParabolicSolver}}} to solve a parabolic time-dependent
 * linear PDE
 *
 * EMPTYLINE
 *
 * The following header files need to be included.
 * First we include the header needed to define this class as a test suite */
#include <cxxtest/TestSuite.h>
/* On some systems there is a clash between Boost Ublas includes and PETSc.  This can be
 * resolved by making sure that Chaste's interface to the Boost libraries are included
 * as early as possible.
 */
#include "UblasCustomFunctions.hpp"
/* This is the class that is needed to solve a linear elliptic pde */
#include "SimpleLinearEllipticSolver.hpp"
/* This is the class that is needed to solve a linear parabolic PDE */
#include "SimpleLinearParabolicSolver.hpp"
/* This is a parabolic PDE, one of the PDEs we will solve */
#include "HeatEquationWithSourceTerm.hpp"
/* We will also solve this PDE */
#include "SimplePoissonEquation.hpp"
/* This is needed to read mesh datafiles of the 'Triangles' format */
#include "TrianglesMeshReader.hpp"
/* This class represents the mesh internally */
#include "TetrahedralMesh.hpp"
/* These are used to specify boundary conditions for the PDEs */
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
/* This class helps us deal with output files */
#include "OutputFileHandler.hpp"
/* !PetscSetupAndFinalize.hpp must be included in every test that uses PETSc. Note that it
 * cannot be included in the source code. */
#include "PetscSetupAndFinalize.hpp"

/* == Test 1: Solving a linear elliptic PDE ==
 *
 * Here, we solve the PDE: div(D grad u) + u + x^2^+y^2^ = 0, in 2D, where
 * D is the diffusion tensor (2 0; 0 1) (ie D11=2, D12=D21=0, D22=1), on a square
 * domain, with boundary conditions u=0 on x=0 or y=0, and (D grad u).n = 0 on x=1 and y=1,
 * where n is the surface normal.
 *
 * EMPTYLINE
 *
 * We need to create a class representing the PDE we want to solve, which will be
 * passed into the solver. The PDE we are solving is of the type
 * {{{AbstractLinearEllipticPde}}}, which is an abstract class with 3 pure methods
 * which have to implemented. The template variables in the following line are both the dimension
 * of the space.
 */
class MyPde : public AbstractLinearEllipticPde<2,2>
{
private:
    /* For efficiency, we will save the diffusion tensor that will be returned by one of the
     * class' methods as a member variable. The diffusion tensor which has to be returned
     * by the {{{GetDiffusionTensor}}} method in PDE classes is of the type
     * {{{c_matrix<double,SIZE,SIZE>}}}, which is a u-blas matrix. We use ublas vectors
     * and matrices where small vectors and matrices are needed. Note that ublas objects
     * are only particularly efficient if optimisation is on ({{{scons build=GccOpt ..}}}).*/
    c_matrix<double,2,2> mDiffusionTensor;

public:
    /* The constructor just sets up the diffusion tensor. We choose a diffusion tensor which
     * corresponds to twice as much diffusion in the x-direction compared to the y-direction */
    MyPde()
    {
        mDiffusionTensor(0,0) = 2.0;
        mDiffusionTensor(0,1) = 0.0;
        mDiffusionTensor(1,0) = 0.0;
        mDiffusionTensor(1,1) = 1.0;
    }

    /* The first method which has to be implemented returns the constant
     * (not dependent on u) part of the source term, which for our PDE is
     * x^2^ + y^2^ */
    double ComputeConstantInUSourceTerm(const ChastePoint<2>& rX)
    {
        return rX[0]*rX[0] + rX[1]*rX[1];
    }

    /* The second method which has to be implemented returns the coefficient in the linear-in-u
     * part of the source term, which for our PDE is just 1.0 */
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>& rX, Element<2,2>* pElement)
    {
        return 1.0;
    }

    /* The third method returns the diffusion tensor D. Note that the diffusion tensor should
     * be symmetric and positive definite for a physical, well-posed problem. */
    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& rX)
    {
        return mDiffusionTensor;
    }
};

/* Next, we define the test suite (a class). It is sensible to name it the same
 * as the filename. The class should inherit from {{{CxxTest::TestSuite}}} */
class TestSolvingLinearPdesTutorial : public CxxTest::TestSuite
{
/* All individual test defined in this test suite '''must''' be declared as public */
public:
    /* Define a particular test */
    void TestSolvingEllipticPde() throw(Exception)
    {
        /* First we declare a mesh reader which reads mesh data files of the 'Triangle'
         * format. The path given is the relative to the main Chaste directory. The reader
         * will look for three datafiles, [name].nodes, [name].ele and (in 2d or 3d)
         * [name].edge. Note that the first template argument here is the dimension of the
         * elements in the mesh ({{{ELEMENT_DIM}}}), and the second is the dimension of the nodes,
         * i.e. the dimension of the space the mesh lives in ({{{SPACE_DIM}}}). Usually
         * {{{ELEMENT_DIM}}} and {{{SPACE_DIM}}} will be equal. */
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        /* Now declare a tetrahedral mesh with the same dimensions */
        TetrahedralMesh<2,2> mesh;
        /* Construct the mesh using the mesh reader */
        mesh.ConstructFromMeshReader(mesh_reader);

        /* Next we instantiate an instance of our PDE we which to solve */
        MyPde pde;

        /* A set of boundary conditions are stored in a {{{BoundaryConditionsContainer}}}. The
         * three template arguments are ELEMENT_DIM, SPACE_DIM and PROBLEM_DIM, the latter being
         * the number of unknowns we are solving for. We have one unknown (ie u is a scalar, not
         * a vector), so in this case {{{PROBLEM_DIM}}}=1. */
        BoundaryConditionsContainer<2,2,1> bcc;

        /* Defining the boundary conditions is the only particularly fiddly part of solving PDEs,
         * unless they are very simple, such as u=0 on the boundary, which could be done
         * as follows */
        //bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        /* We want to specify u=0 on x=0 and y=0. To do this, get a boundary node iterator
         * from the mesh */
        TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        /* Then loop over the boundary nodes, getting the x and y value */
        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            /* if x=0 or y=0.. */
            if ((x==0) || (y==0))
            {
                /* ..create a new {{{ConstBoundaryConditions}}} object. This is a subclass of
                 * {{{AbstractBoundaryCondition}}}, and tells the caller what value to return
                 * given a particular point in space. In the first line below we say that value
                 * should be 0.0. The second line tells the {{{BoundaryConditionsContainer}}}
                 * object that it should associate this boundary condition with this node
                 * ({{{*iter}}} being a pointer to a {{{Node<2>}}}).
                 */
                ConstBoundaryCondition<2>* p_dirichlet_boundary_condition = new ConstBoundaryCondition<2>(0.0);
                bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            }
            iter++;
        }

        /* Now we create Neumann boundary conditions for the ''surface elements'' on x=1 and y=1. Note that
         * Dirichlet boundary conditions are defined on nodes, whereas Neumann boundary conditions are
         * defined on surface elements. Note also that the natural boundary condition statement for this
         * PDE is (D grad u).n = g(x) (where n is the outward-facing surface normal), and g(x) is a prescribed
         * function, ''not'' something like du/dn=g(x). Hence the boundary condition we are specifying is
         * (D grad u).n = 0.
         *
         * EMPTYLINE
         *
         * '''Important note for 1D:''' This means that if we were solving 2u,,xx,,=f(x) in 1D, and
         * wanted to specify du/dx=1 on the LHS boundary, the Neumann boundary value we have to specify is
         * -2, as n=-1 (outward facing normal) so (D gradu).n = -2 when du/dx=1.
         *
         * EMPTYLINE
         *
         * To define Neumann bcs, we define another constant boundary condition object (created using
         * {{{new}}} - note that the {{{BoundaryConditionsContainer}}} object deals with deleting
         * its {{{AbstractBoundaryCondition}}} objects), and then loop over surface elements, using the
         * iterator provided by the mesh class.
         */
        ConstBoundaryCondition<2>* p_neumann_boundary_condition = new ConstBoundaryCondition<2>(0.0);

        TetrahedralMesh<2,2>::BoundaryElementIterator surf_iter
          = mesh.GetBoundaryElementIteratorBegin();
        while (surf_iter < mesh.GetBoundaryElementIteratorEnd())
        {
            /* Get the x and y values of any node (here, the 0th) in the element */
            unsigned node_index = (*surf_iter)->GetNodeGlobalIndex(0);
            double x = mesh.GetNode(node_index)->GetPoint()[0];
            double y = mesh.GetNode(node_index)->GetPoint()[1];

            /* if x=1 or y=1.. */
            if ( (fabs(x-1.0) < 1e-6) || (fabs(y-1.0) < 1e-6) )
            {
                /* associate the boundary condition with the surface element */
                bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
            }

            /* and increment the iterator */
            surf_iter++;
        }

        /* Next we define the solver of the PDE.
         * To solve {{{AbstractLinearEllipticPde}}} (which is the type of pde {{{MyPde}}} is),
         * we use a {{{SimpleLinearEllipticSolver}}}. The solver, again templated over
         * {{{ELEMENT_DIM}}} and {{{SPACE_DIM}}}, needs to be given (pointers to) the mesh,
         * pde and boundary conditions.
         */
        SimpleLinearEllipticSolver<2,2> solver(&mesh,&pde,&bcc);

        /* To solve, just call {{{Solve()}}}. A Petsc vector is returned. */
        Vec result = solver.Solve();

        /* It is a pain to access the individual components of a Petsc vector, even in
         * sequential. A helper class called {{{ReplicatableVector}}} has been created. Create
         * an instance of one of these, using the Petsc {{{Vec}}} as the data. The ith
         * component of {{{result}}} can now be obtained by simply doing {{{result_repl[i]}}}.
         */
        ReplicatableVector result_repl(result);

        /* Let us write out the solution to a file. To do this, create an
         * {{{OutputFileHandler}}}, passing in the directory we want files written to.
         * This is relative to the directory defined by the CHASTE_TEST_OUTPUT environment
         * variable - usually /tmp/chaste/testoutput. Note by default the output directory
         * passed in is cleaned. To avoid this, {{{false}}} can be passed in as a second
         * parameter
         */
        OutputFileHandler output_file_handler("TestSolvingLinearPdeTutorial");

        /* Create an {{{out_stream}}}, which is a stream to a particular file. An {{{out_stream}}}
         * is a pointer to a ofstream */
        out_stream p_file = output_file_handler.OpenOutputFile("linear_solution.txt");

        /* Loop over the entries of the solution */
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            /* Get the x and y-values of the node corresponding to this entry. The method
             * {{{GetNode}}} on the mesh class returns a pointer to a {{{Node}}} */
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            /* Get the computed solution at this node from the {{{ReplicatableVector}}} */
            double u = result_repl[i];

            /* Finally, write x, y and u to the output file. The solution could then be
             * visualised in (eg) matlab, using the commands:
             * {{{sol=load('linear_solution.txt'); plot3(sol(:,1),sol(:,2),sol(:,3),'.');}}}*/
            (*p_file) << x << " " << y << " " << u << "\n";
        }

        /* All Petsc {{{Vec}}}s should be destroyed when they are no longer needed */
        VecDestroy(result);
    }

    /*
     * == Test 2: Solving a linear parabolic PDE ==
     *
     * Now we solve a parabolic PDE. We choose a simple problem so that the code changes
     * needed from the elliptic case are clearer. We will solve
     * du/dt = div(grad u) + u, in 3d, with boundary conditions u=1 on the boundary, and initial
     * conditions u=1
     *
     */
    void TestSolvingParabolicPde() throw(Exception)
    {
        /* Create a 10 by 10 by 10 mesh in 3D, this time using the {{{ConstructRegularSlabMesh}}} method
         * on the mesh. The first parameter is the cartesian space-step and the other three parameters are the width, height and depth of the mesh.*/
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructRegularSlabMesh(0.1, 1.0, 1.0, 1.0);

        /* Our PDE object should be a class that is derived from the {{{AbstractLinearParabolicPde}}}.
         * We could write it ourselves as in the previous test, but since the PDE we want to solve is
         * so simple, it has already been defined (look it up! - it is located in pde/test/pdes)
         */
        HeatEquationWithSourceTerm<3> pde;

        /* Create a new boundary conditions container and specify u=1.0 on the boundary */
        BoundaryConditionsContainer<3,3,1> bcc;
        bcc.DefineConstantDirichletOnMeshBoundary(&mesh, 1.0);

        /* Create an instance of the solver, passing in the mesh, pde and boundary conditions.
         */
        SimpleLinearParabolicSolver<3,3> solver(&mesh,&pde,&bcc);

        /* For parabolic problems, initial conditions are also needed. The solver will expect
         * a Petsc vector, where the i-th entry is the initial solution at node i, to be passed
         * in. To create this Petsc {{{Vec}}}, we will use a helper function in the {{{PetscTools}}}
         * class to create a {{{Vec}}} of size num_nodes, with each entry set to 1.0. Then we
         * set the initial condition on the solver */
        Vec initial_condition = PetscTools::CreateAndSetVec(mesh.GetNumNodes(), 1.0);
        solver.SetInitialCondition(initial_condition);

        /* Next define the start time, end time, and timestep, and set them. */
        double t_start = 0;
        double t_end = 1;
        double dt = 0.01;
        solver.SetTimes(t_start, t_end);
        solver.SetTimeStep(dt);

        /* Now we can solve the problem. The {{{Vec}}} that is returned can be passed into a
         * {{{ReplicatableVector}}} as before
         */
        Vec solution = solver.Solve();
        ReplicatableVector solution_repl(solution);

        /* Let's also solve the equivalent static PDE, ie set du/dt=0, so 0=div(gradu) + u. This
         * is easy, as the PDE class has already been defined */
        SimplePoissonEquation<3,3> static_pde;
        SimpleLinearEllipticSolver<3,3> static_solver(&mesh, &static_pde, &bcc);
        Vec static_solution = static_solver.Solve();
        ReplicatableVector static_solution_repl(static_solution);

        /* We can now compare the solution of the parabolic PDE at t=1 with the static solution,
         * to see if the static equilibrium solution was reached in the former. (Ideally we should
         * compute some relative error, but we just compute an absolute error for simplicity). */
        for (unsigned i=0; i<static_solution_repl.GetSize(); i++)
        {
            TS_ASSERT_DELTA( solution_repl[i], static_solution_repl[i], 1e-3);
        }

        /* All Petsc vectors should be destroyed when they are no longer needed */
        VecDestroy(initial_condition);
        VecDestroy(solution);
        VecDestroy(static_solution);
    }
};

#endif /*TESTSOLVINGLINEARPDESTUTORIAL_HPP_*/
