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
 *
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 *
 *
 *
 */
#ifndef TESTWRITINGPDESOLVERSTUTORIAL_HPP_
#define TESTWRITINGPDESOLVERSTUTORIAL_HPP_

#include <cxxtest/TestSuite.h>



/*
 * == Introduction ==
 *
 * Chaste can be used to set up solvers for more general (coupled) PDEs, for which the
 * user just, essentially, needs to code up the integrands of any finite element (FE) matrices
 * or vectors, without having to deal with set-up, looping over elements, numerical
 * quadrature, assembly, or solving the linear system. If you have a coupled set of
 * linear PDEs for which it is appropriate to use linear basis functions for each unknown
 * (for example, a coupled set of reaction-diffusion equations), then it is relatively
 * straightforward to set up a solver which will be parallel and reliable (as all the base
 * components are heavily tested).
 *
 * EMPTYLINE
 *
 * There are solvers for general simple (uncoupled) linear PDEs already provided, such
 * as the `SimpleLinearEllipticSolver`. These are for PDEs that can be written in a generic
 * form (`SimpleLinearEllipticPde`, for example). However more general (coupled) PDEs can't be
 * written a generic form, so the user has to write their own solver.
 *
 * EMPTYLINE
 *
 * For this tutorial the user certainly ought to have read the solving-PDEs
 * tutorials. Also, it is helpful to read the associated lectures nodes: maybe the
 * slides on solving equations using finite elements if you are not familiar with this
 * (lecture 2), but especially the slides on the general design of the Chaste finite element
 * solvers (lecture 3).
 *
 * EMPTYLINE
 *
 * Let us use the terminology "assembled in an FE manner" for any matrix or vector that
 * defined via a volume/surface/line integral, and which is constructed by: looping over
 * elements (or surface elements, etc), computing the elemental contribution (ie a small
 * matrix/vector) using numerical quadrature, and adding to the full matrix/vector.
 *
 * EMPTYLINE
 *
 * We only consider linear problems here. In these problems the discretised FE problem leads to
 * a linear system, Ax=b, to be solved (at each timestep in time-dependent problems).
 * There are two cases to be distinguished. The first case is where BOTH A and b are
 * 'assembled in an FE manner', b possibly being composed of a volume integral plus a
 * surface integral. The other case is where this is not true, for example where
 * b = Mc+d, where d (vector) and M (matrix) are assembled in an FE manner, but not c.
 *
 * EMPTYLINE
 *
 * The chaste PDE classes include ASSEMBLER classes (for setting up anything assembled in an FE manner),
 * and SOLVER classes (for setting up linear systems). In the general case, solvers need to own assemblers
 * for setting up each part of the linear system. However for the former case (both A and b in Ax=b are
 * assembled), we can use the design where the solver IS AN assembler. We illustrate how to do
 * this in this first tutorial
 *
 * == Writing solvers ==
 *
 * Let us write a solver for the coupled 2-unknown problem
 * {{{
 *    Laplacian(u) + v = f(x,y)
 *    Laplacian(v) + u = g(x,y)
 * }}}
 * (`Laplacian(u)` of course representing u_xx_+u_yy_), and
 * where f and g are chosen so that (with zero-dirichlet boundary conditions)
 * the solution is: u = sin(pi*x)sin(pi*x), v = sin(2*pi*x)sin(2*pi*x)
 *
 * EMPTYLINE
 *
 * Using linear basis functions, and a mesh with N nodes, the linear system that needs to be set up is
 * of size 2N by 2N, and in block form is:
 * {{{
 *   [ K   -M  ] [U]  =  [b1]
 *   [ -M   K  ] [V]     [b2]
 * }}}
 * where `K` is the stiffness matrix, `M` the mass matrix, `U` the vector of nodal values
 * of u, `V` the vector of nodal values of v, `b1_i = integral(f\phi_i dV)` and
 * `b1_i = integral(g\phi_i dV)`, where the basis functions are `phi_i`.
 *
 * This is the linear system which we now write a solver to set up. Note, however, that
 * the main Chaste solvers assume a STRIPED data format, ie that the unknown vector
 * is `[U_1 V_1 U_2 V_2 .. U_n V_n]`, not `[ U_1 U_2 .. U_n V_1 V_2 .. V_n]`. We write down
 * equations in block form as it makes things clearer, but have to remember that the code
 * deals with striped data structures.
 *
 * EMPTYLINE
 *
 * These are some basic includes as in the solving-PDEs tutorial
 */
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
/* These two classes will be used by the solver */
#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
/* We will solve a second problem, below, which will be time-dependent and will
 * require this */
#include "AbstractDynamicLinearPdeSolver.hpp"

/* The linear system is Ax=b where A and b are FE assembled, so we can use the solver-is-an-assembler
 * design. To construct our solver, we inherit from `AbstractAssemblerSolverHybrid` which links
 * the solver clases to the assembler classes, and `AbstractStaticLinearPdeSolver`. (For time-dependent
 * problems, the second parent would be `AbstractDynamicLinearPdeSolver`). Note the template
 * parameter `PROBLEM_DIM` below, in this case it is 2 as there are two unknowns.
 */
class MyCoupledStaticPdeSolver
    : public AbstractAssemblerSolverHybrid<2/*elem_dim*/,2/*space_dim*/,2/*problem_dim*/,NORMAL/*amount of interpolation*/>,
      public AbstractStaticLinearPdeSolver<2,2,2>
{
private:
    /* The function f */
    double f(double x,double y)
    {
        return -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y) + sin(2*M_PI*x)*sin(2*M_PI*y);
    }
    /* The function g */
    double g(double x,double y)
    {
        return -8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y) + sin(M_PI*x)*sin(M_PI*y);
    }

    /*
     *  The abstract assembler parent classes know how to assemble matrices and vectors, but the concrete
     *  class needs to provide the integrand of the elemental contribution to A and b. This first
     *  method returns the elemental contribution, given the provided bases (`rPhi`, `rGradPhi`)
     *  of the matrix A. The '3's here represent the number of bases per element (ie the number
     *  of nodes as linear bases are being used).
     */
    c_matrix<double,2*3,2*3> ComputeMatrixTerm(c_vector<double, 3>& rPhi /* the three bases for the current element, evaluated at the current quad pt*/,
                                               c_matrix<double, 2, 3>& rGradPhi /* gradients of the three bases */,
                                               ChastePoint<2>& rX           /* phsyical coordinate of quad point */,
                                               c_vector<double,2>& rU       /* current solution (unused here as a linear static problem */,
                                               c_matrix<double,2,2>& rGradU /* current solution gradient (unused here as a linear static problem */,
                                               Element<2,2>* pElement)
    {
        /*
         * Set up the matrix, which corresponds to the elemental contribution for the matrix
         * written above, taking into account the striped nature of the matrices and vectors.
         * (Note: the following can be done more efficiently using matrix slices and products,
         * see `BidomainAssembler` for example).
         */
        c_matrix<double,2*(2+1),2*(2+1)> ret = zero_matrix<double>(2*(2+1), 2*(2+1));

        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                for (unsigned k=0; k<2; k++)
                {
                    ret(2*i,  2*j)   += rGradPhi(k,i)*rGradPhi(k,j);
                    ret(2*i+1,2*j+1) += rGradPhi(k,i)*rGradPhi(k,j);
                }

                ret(2*i+1, 2*j)   = -rPhi(i)*rPhi(j);
                ret(2*i,   2*j+1) = -rPhi(i)*rPhi(j);
            }
        }
        return ret;
    }

    /* Similarly compute the elemental contribution to the RHS vector */
    c_vector<double,2*(2+1)> ComputeVectorTerm(c_vector<double, 2+1>& rPhi,
                                               c_matrix<double, 2, 2+1>& rGradPhi,
                                               ChastePoint<2>& rX,
                                               c_vector<double,2>& rU,
                                               c_matrix<double,2,2>& rGradU,
                                               Element<2,2>* pElement)
    {
       c_vector<double,2*(2+1)> ret;

        for (unsigned i=0; i<3; i++)
        {
            ret(2*i)   = -f(rX[0],rX[1]) * rPhi(i);
            ret(2*i+1) = -g(rX[0],rX[1]) * rPhi(i);
        }
        return ret;
    }
    /* NOTE: we will not be solving this equation subject to any non-zero Neumann
     * boundary conditions. If we were though, we would have to also provide a method
     * `ComputeSurfaceVectorTerm(..)`.
     *
     * EMPTYLINE
     *
     * These classes which inherit from both assemblers and solvers must
     * provide the following method, which links the two. Just copy and paste.
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
    }

public:
    /* The constructor takes in a mesh and boundary conditions container, and passes
     * them to the parent classes.
     */
    MyCoupledStaticPdeSolver(TetrahedralMesh<2,2>* pMesh,
                             BoundaryConditionsContainer<2,2,2>* pBoundaryConditions)
        : AbstractAssemblerSolverHybrid<2,2,2,NORMAL>(pMesh,pBoundaryConditions),
          AbstractStaticLinearPdeSolver<2,2,2>(pMesh)
    {
    }
};


/*
 *  That is the solver written. To solve it, we do the following, which is more-or-less
 *  the same as with the other PDE solvers from the previous tutorials.
 *
 */
class TestWritingPdeSolversTutorial : public CxxTest::TestSuite
{
public:
    void TestRealCoupledPde() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,0); // zero dirichlet for u
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,1); // zero dirichlet for v

        // Purpose-made solver for this problem:
        MyCoupledStaticPdeSolver solver(&mesh,&bcc);

        /* The `AbstractStaticLinearPdeSolver` class from which our solver
         * inherits, provides a `Solve` method.
         */
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        /* Compare against the exact solution
         */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];

            double u = result_repl[2*i];
            double v = result_repl[2*i+1];

            double u_exact = sin(M_PI*x)*sin(M_PI*y);
            double v_exact = sin(2*M_PI*x)*sin(2*M_PI*y);

            TS_ASSERT_DELTA( u, u_exact, 0.002);
            TS_ASSERT_DELTA( v, v_exact, 0.007);
        }

        VecDestroy(result);
    }
};

#endif // TESTWRITINGPDESOLVERSTUTORIAL_HPP_
