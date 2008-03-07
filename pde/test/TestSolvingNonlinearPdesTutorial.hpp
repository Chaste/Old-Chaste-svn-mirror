/*
 * 
 * 
 * 
 *  Chaste tutorial - this page get automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 * 
 * 
 * 
 * 
 * 
 */  
#ifndef TESTSOLVINGNONLINEARPDESTUTORIAL_HPP_
#define TESTSOLVINGNONLINEARPDESTUTORIAL_HPP_

/* 
 * In this tutorial we show how Chaste can be used to solve a nonlinear elliptic PDEs. 
 * We will solve the PDE div.(u grad u) + 1 = 0, on a square domain, with boundary
 * conditions u=0 on y=0; and Neumann boundary conditions: (u grad u).n = 0 on x=0 and x=1;
 * and (u grad u).n = y on y=1. 
 * 
 * EMPTYLINE
 * 
 * For nonlinear PDEs, the finite element equations are of the form F(U)=0, where
 * U=(U,,1,, , ... , U,,N,,) is a vector of the unknowns at each node, and F some 
 * non-linear vector valued function. To solve this, a nonlinear solver is required. 
 * Chaste can solve this with Newton's method, or (default) Petsc's nonlinear solvers. 
 * Solvers of such nonlinear problems usually require the Jacobian of the problem, ie the
 * matrix A = dF/dU, or at least an approximation of the Jacobian. 
 * 
 * EMPTYLINE
 * 
 * The following header files need to be included, as in the linear PDEs tutorial
 */
#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "TrianglesMeshReader.cpp"
#include "PetscSetupAndFinalize.hpp"
/* This is the assembler for nonlinear elliptic PDEs */
#include "SimpleNonlinearEllipticAssembler.hpp"
/* In this test we also show how to define Neumman boundary conditions which 
 * depend on spatial location, for which the following class is needed */
#include "FunctionalBoundaryCondition.hpp"
/* We will choose to use the Chaste newton solver rather than Petsc's nonlinear
 * solver */ 
#include "SimpleNewtonNonlinearSolver.hpp"

/* As in the linear PDEs tutorial, we have to define the PDE class we want to
 * solve (assuming one has not already been created). Nonlinear elliptic PDEs
 * should inherit from {{{AbstractNonlinearEllipticPde}}}, which has five pure
 * methods which have to be implemented in this concrete class. Here, we define 
 * the PDE div.(u grad u) + 1 = 0.
 */
class MyNonlinearPde : public AbstractNonlinearEllipticPde<2>
{
public:
    /* The first is the part of the source term that is independent of u */
    double ComputeLinearSourceTerm(const ChastePoint<2>& rX)
    {
        return 1.0;
    }
    
    /* The first is the part of the source term that is dependent on u */
    double ComputeNonlinearSourceTerm(const ChastePoint<2>& rX, double u)
    {
        return 0.0;
    }
    
    /* The third is the diffusion tensor, which unlike the linear case, can be
     * dependent on u */
    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& rX, double u)
    {
        return identity_matrix<double>(2)*u;
    }
    
    /* We also need to provide the derivatives with respect to u of the last two methods,
     * so that the Jacobian matrix can be assembled. The derivative of the nonlinear source
     * term is
     */ 
    double ComputeNonlinearSourceTermPrime(const ChastePoint<2>& , double )
    {
        return 0.0;
    }

    /* And the derivative of the diffusion tensor is just the identity matrix */   
    c_matrix<double,2,2> ComputeDiffusionTermPrime(const ChastePoint<2>& rX, double u)
    {
        return identity_matrix<double>(2);
    }
};

/* We also need to define a (global) function that we become the Neumman boundary 
 * conditions, via the {{{FunctionalBoundaryCondition}}} class (see below). This
 * function is f(x,y) = y
 */
double MyNeummanFunction(const ChastePoint<2>& rX)
{
    return rX[1];
}
 

/* Next, we define the test suite, as before */
class TestSolvingNonlinearPdesTutorial : public CxxTest::TestSuite
{
public:
    /* Define a particular test. Note the {{{throw(Exception)}}} at the end of the
     * declaration. This causes {{{Exception}}} messages to be printed out if an
     * {{{Exception}}} is thrown, rather than just getting the message "terminate 
     * called after throwing an instance of 'Exception' " */
    void TestSolvingNonlinearEllipticPde() throw(Exception)
    {
        /* As usual, first create a mesh */
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        /* Next, instantiate the PDE to be solved */
        MyNonlinearPde pde;
        
        /* 
         * Then we have to define the boundary conditions. First, the Dirichlet boundary
         * condition, u=0 on x=0, using the boundary node iterator
         */
        BoundaryConditionsContainer<2,2,1> bcc;         
        ConstBoundaryCondition<2>* p_zero_bc = new ConstBoundaryCondition<2>(0.0);
        for( ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
             node_iter != mesh.GetBoundaryNodeIteratorEnd();
             node_iter++)
        {
            if (fabs((*node_iter)->GetPoint()[1]) < 1e-12)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, p_zero_bc);
            }
        }
        
        /* And then the Neumman conditions. Neumann boundary condition are defined on
         * surface elements, and for this problem, the Neumman boundary value depends
         * on the position in space, so we make use of the {{{FunctionalBoundaryCondition}}}
         * object, which contains a pointer to a function, and just returns the value 
         * of that function for the required point when the {{{GetValue}}} method is called.
         */ 
        FunctionalBoundaryCondition<2>* p_functional_bc
          = new FunctionalBoundaryCondition<2>( &MyNeummanFunction );
        /* Next, loop over surface elements */
        for( ConformingTetrahedralMesh<2,2>::BoundaryElementIterator elt_iter = mesh.GetBoundaryElementIteratorBegin();
             elt_iter != mesh.GetBoundaryElementIteratorEnd();
             elt_iter++ )
        {
            /* Get the y value of any node (here, the zero-th) */
            double y = (*elt_iter)->GetNodeLocation(0,1);
            /* If y=1.. */
            if (fabs(y-1.0) < 1e-12)
            {
                /* .. then associate the functional boundary condition, (Dgradu).n = y,
                 *  with the surface element */
                bcc.AddNeumannBoundaryCondition(*elt_iter, p_functional_bc);
            }
            else
            {
                /* else associate the zero boundary condition (ie zero flux) with this
                 * element */
                bcc.AddNeumannBoundaryCondition(*elt_iter, p_zero_bc);
            }
        }
        /* Note that in the above loop, the zero Neumman boundary condition was applied
         * to all surface elements for which y!=1, which included the Dirichlet surface
         * y=0. This is OK, as Dirichlet boundary conditions are applied to the finite
         * element matrix after Neumman boundary conditions, where the appropriate rows
         * in the matrix are overwritten
         * 
         * EMPTYLINE 
         * 
         * This is the assembler for solving nonlinear problems, which, as usual,
         * takes in the mesh, the PDE, and the boundary conditions */
        SimpleNonlinearEllipticAssembler<2,2> assembler(&mesh, &pde, &bcc);
        
        /* The assembler also needs to be given an initial guess, which will be
         * a Petsc vector. We can make use of a helper method to create it. (We
         * could also have used the {{{PetscTools}}} class as in the previous
         * tutorial) */ 
        Vec initial_guess = assembler.CreateConstantInitialGuess(0.25);
        
        /* '''Optional:''' To use Chaste's Newton solver to solve nonlinear vector equations that are
         * assembled, rather than the default Petsc nonlinear solvers, we can
         * do the following: */
        SimpleNewtonNonlinearSolver newton_solver;
        assembler.SetNonlinearSolver(&newton_solver);
        /* '''Optional:''' We can also manually set tolerances, and whether to print statistics, with
         * this nonlinear vector equation solver */
        newton_solver.SetTolerance(1e-10);
        newton_solver.SetWriteStats();

        /* Now call {{{Solve}}}, passing in the initial guess */
        Vec answer = assembler.Solve(initial_guess);

        /* Note that we could have got the assembler to not use an analytical Jacobian
         * and use a numerically-calculated Jacobian instead, by passing in false as a second 
         * parameter
         */
        //Vec answer = assembler.Solve(initial_guess, false);
        
        
        /* Once solved, we can check the obtained solution against the analytical 
         * solution */
        ReplicatableVector answer_repl(answer);
        for (unsigned i=0; i<answer_repl.size(); i++)
        {
            double y = mesh.GetNode(i)->GetPoint()[1];
            double exact_u = sqrt(y*(4-y));
            TS_ASSERT_DELTA(answer_repl[i], exact_u, 0.15);
        }

        /* Finally, we have to remember to destroy the Petsc {{{Vec}}}s */
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
};

#endif /*TESTSOLVINGNONLINEARPDESTUTORIAL_HPP_*/
