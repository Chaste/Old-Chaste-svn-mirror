#ifndef _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_



#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include <vector>
#include <cmath>
//#include <iostream>

#include "SimpleNonlinearEllipticAssembler.hpp"
#include "SimplePetscNonlinearSolver.hpp"

#include "BoundaryConditionsContainer.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "ConstBoundaryCondition.hpp"

#include "NonlinearHeatEquationPde.hpp"
#include "NonlinearHeatEquation2Pde.hpp"
#include "NonlinearHeatEquation3Pde.hpp"
#include "NonlinearHeatEquation4Pde.hpp"
#include "NonlinearHeatEquation5Pde.hpp"
#include "Example2DNonlinearEllipticPde.hpp"
#include "NonlinearLinearHeatEquationPde.hpp"
#include "ExampleNasty2dNonlinearEllipticPde.hpp"
#include "TrianglesMeshReader.cpp"
#include "PetscSetupAndFinalize.hpp"



/**
 * For use in TestSimpleNonlinearEllipticAssembler::Test2dOnUnitSquare.
 */
double bc_x1_func(Point<2> p)
{
    return 2*(2+p[1]*p[1]);
}
/**
 * For use in TestSimpleNonlinearEllipticAssembler::Test2dOnUnitSquare.
 */
double bc_y1_func(Point<2> p)
{
    return 2*(2+p[0]*p[0]);
}

/**
 * For use in TestSimpleNonlinearEllipticAssembler::TestNasty2dEquationOnUnitSquare.
 */
double bc_x1_func2(Point<2> p)
{
    return sin(2)*(sin(1)*sin(1)+1+p[1]*p[1]);
}
/**
 * For use in TestSimpleNonlinearEllipticAssembler::TestNasty2dEquationOnUnitSquare.
 */
double bc_y1_func2(Point<2> p)
{
    return 2*(2+sin(p[0])*sin(p[0]));
}

/**
 * For use in TestSimpleNonlinearEllipticAssembler::TestWithHeatEquation2DAndNeumannBCs
 */
double one_bc(Point<2> p)
{
    return p[1];
}


class TestSimpleNonlinearEllipticAssembler : public CxxTest::TestSuite
{

private:

    /**
     * Refactor code to set up a PETSc vector holding the initial guess.
     */
    Vec CreateInitialGuessVec(int size)
    {
        Vec initial_guess;
        VecCreate(PETSC_COMM_WORLD, &initial_guess);
        VecSetSizes(initial_guess, PETSC_DECIDE, size);
        VecSetFromOptions(initial_guess);
        return initial_guess;
    }
   
    
public:

    void TestAssembleResidual( void )
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/practical1_1d_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        //Adding Dirichlet BC at node 0
        double DirichletBCValue = 5.0;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(DirichletBCValue);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        
        // adding von Neumann BC at the last node
        double VonNeumannBCValue = 9.0;
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(VonNeumannBCValue);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--; // to be consistent with c++ :))), GetBoundaryElementIteratorEnd points to one element passed it
        bcc.AddNeumannBoundaryCondition(*iter,pBoundaryCondition1);
        
        // initialize 'solution' vector
        Vec solution;
        VecCreate(PETSC_COMM_WORLD, &solution);
        VecSetSizes(solution, PETSC_DECIDE, mesh.GetNumNodes());
        VecSetFromOptions(solution);
        
        NonlinearHeatEquationPde<1> pde;
        
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);
                
        // Set 'solution' to 1 and compute residual
        double h = 0.01;
        for (int global_index = 0; global_index<mesh.GetNumNodes(); global_index++)
        {
            VecSetValue(solution, global_index, (PetscReal) 1, INSERT_VALUES);
        }
        double initial_guess = 1.0;
        VecAssemblyBegin(solution);
        VecAssemblyEnd(solution);
        Vec residual;
        VecDuplicate(solution, &residual);
        
        assembler.AssembleResidual(solution, residual);
        
        PetscScalar *p_residual;
        VecGetArray(residual, &p_residual);
        
        int lo, hi;
        VecGetOwnershipRange(residual, &lo, &hi);
        
        if (lo<=0 && 0<hi)
        {
            int local_index = 0-lo;
            double value1 = p_residual[local_index];
            TS_ASSERT(fabs(value1 + DirichletBCValue - initial_guess) < 0.001);
        }
        if (lo<=1 && 1<hi)
        {
            int local_index = 1-lo;
            double value2 = p_residual[local_index];
            TS_ASSERT(fabs(value2 + h) < 0.001);
        }
        if (lo<=mesh.GetNumNodes()-1 && mesh.GetNumNodes()-1<hi)
        {
            int local_index = mesh.GetNumNodes()-1-lo;
            double valueLast = p_residual[local_index];
            TS_ASSERT(fabs(valueLast + VonNeumannBCValue + h/2) < 0.001);
        }
        VecRestoreArray(residual, &p_residual);
        VecDestroy(residual);
        VecDestroy(solution);
    }
    
    void TestNumericalAgainstAnalyticJacobian()
    {
        PetscInt n = 11;  // Mesh size
        Mat numerical_jacobian;
#if (PETSC_VERSION_MINOR == 2) //Old API
        MatCreate(PETSC_COMM_WORLD, PETSC_DETERMINE, PETSC_DETERMINE, n, n, &numerical_jacobian);
#else
        MatCreate(PETSC_COMM_WORLD, &numerical_jacobian);
        MatSetSizes(numerical_jacobian, PETSC_DETERMINE, PETSC_DETERMINE, n, n);
#endif
        //MatSetType(numerical_jacobian, MATSEQDENSE);
        MatSetFromOptions(numerical_jacobian);
        Mat analytic_jacobian;
#if (PETSC_VERSION_MINOR == 2) //Old API
        MatCreate(PETSC_COMM_WORLD, PETSC_DETERMINE, PETSC_DETERMINE, n, n, &analytic_jacobian);
#else
        MatCreate(PETSC_COMM_WORLD, &analytic_jacobian);
        MatSetSizes(analytic_jacobian, PETSC_DETERMINE, PETSC_DETERMINE, n, n);
#endif
        //MatSetType(analytic_jacobian, MATSEQDENSE);
        MatSetFromOptions(analytic_jacobian);
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        NonlinearHeatEquationPde<1> pde;


        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), p_boundary_condition);


        // assembler
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);


        // Set up initial solution guess for residuals
        int length=mesh.GetNumNodes();
        Vec initial_guess;
        VecCreate(PETSC_COMM_WORLD, &initial_guess);
        VecSetSizes(initial_guess, PETSC_DECIDE, length);
        VecSetFromOptions(initial_guess);
        for (int i=0; i<length ; i++)
        {
            VecSetValue(initial_guess, i, (-0.01*i*i), INSERT_VALUES);
        }
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);

        int errcode = assembler.AssembleJacobianNumerically(initial_guess, &numerical_jacobian);
        TS_ASSERT_EQUALS(errcode, 0);
        
        assembler.mUseAnalyticalJacobian = true; // can access the member variable as this class is a friend
        errcode = assembler.AssembleJacobian(initial_guess, &analytic_jacobian);
        TS_ASSERT_EQUALS(errcode, 0);
        MatAssemblyBegin(numerical_jacobian, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(numerical_jacobian, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(analytic_jacobian, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(analytic_jacobian, MAT_FINAL_ASSEMBLY);
//		TS_TRACE("Numerical:");
//		MatView(numerical_jacobian, 0);
//		TS_TRACE("Analytical:");
//		MatView(analytic_jacobian, 0);

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
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearHeatEquationPde<1> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), p_boundary_condition);
        
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);

        // Set up initial guess
        Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            VecSetValue(initial_guess, i, (-0.01*i*i), INSERT_VALUES);
        }
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);
        
        Vec answer = assembler.Solve(initial_guess, true);
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double u = sqrt(x*(1-x));
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.001);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void TestWithHeatEquation1DAndNeumannBCs()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearHeatEquationPde<1> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        // u(0) = 0
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), p_boundary_condition);
        // u(1)*u'(1) = 1
        p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
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
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double u = sqrt(x*(4-x));
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.001);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void TestWithHeatEquation1D2()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearHeatEquation2Pde<1> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), p_boundary_condition);
        p_boundary_condition = new ConstBoundaryCondition<1>(exp(1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), p_boundary_condition);
        
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);
        
        // Set up initial Guess
        Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            VecSetValue(initial_guess, i, (1.0+0.01*i*i), INSERT_VALUES);
        }
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);
        
        Vec answer = assembler.Solve(initial_guess, true);
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double u = exp(0.5*(3.0*x-x*x));
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.001);
        }
        
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void TestWithHeatEquation1D3()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearHeatEquation3Pde<1> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(sqrt(2.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), p_boundary_condition);
        p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), p_boundary_condition);
        
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc,3);
        
        // Set up initial Guess
        Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
        for (int global_index=0; global_index<mesh.GetNumNodes(); global_index++)
        {
            VecSetValue(initial_guess, global_index, (1.5-0.15*global_index), INSERT_VALUES);
        }
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);
        
        Vec answer = assembler.Solve(initial_guess, true);
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double u = sqrt(2.0*(exp(-x)-x*exp(-1.0)));
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.001);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void TestWithHeatEquation1D4()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearHeatEquation4Pde<1> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        // u(1) = exp(1.0)
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(exp(-1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), p_boundary_condition);
        // u(0)^2*u'(0) = 0.0
        p_boundary_condition = new ConstBoundaryCondition<1>(0.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
        
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);
        
        // Set up initial Guess
        Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
        double x1;
        for (int global_index=0; global_index<mesh.GetNumNodes(); global_index++)
        {
            x1=0.1*(double)(global_index);
            VecSetValue(initial_guess, global_index, 0.35*(1-x1*x1), INSERT_VALUES);
        }
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);
        
        Vec answer = assembler.Solve(initial_guess, true);
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double u = x*exp(-x);
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.01);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void TestWithHeatEquation1D5()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearHeatEquation5Pde<1> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        // u(1) = exp(-1.0)
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(exp(-1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), p_boundary_condition);
        // u(0)^2*u'(0) = -1.0
        // Note that we specify 1 as the value, since figuring out which direction
        // the normal is in is hard in 1D.
        p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
        
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);
        
        // Set up initial Guess
        Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
        double x1;
        for (int global_index=0; global_index<mesh.GetNumNodes(); global_index++)
        {
            x1=0.1*(double)(global_index);
            VecSetValue(initial_guess, global_index, 0.35*(1-x1*x1), INSERT_VALUES);
        }
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);
        
        Vec answer = assembler.Solve(initial_guess, true);
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double u = exp(-x);
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.01);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void TestWithHeatEquation1DAndNeumannBCs2()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearHeatEquationPde<1> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        // u(1) = sqrt(3)
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(sqrt(3));
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), p_boundary_condition);
        // u(0)*u'(0) = 2
        // Note that we specify -2 as the value, since figuring out which direction
        // the normal is in is hard in 1D.
        p_boundary_condition = new ConstBoundaryCondition<1>(-2.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
        
        SimpleNonlinearEllipticAssembler<1,1> assembler(&mesh, &pde, &bcc);
        
        // cover the bad size exception
        Vec badly_sized_init_guess = CreateInitialGuessVec(1);
        double value = 1;
#if (PETSC_VERSION_MINOR == 2) //Old API
        VecSet(&value, badly_sized_init_guess);
#else
        VecSet(badly_sized_init_guess, value);
#endif
        VecAssemblyBegin(badly_sized_init_guess);
        VecAssemblyEnd(badly_sized_init_guess);
        
        TS_ASSERT_THROWS_ANYTHING( assembler.Solve(badly_sized_init_guess, true) );
        
        // Set up initial Guess
        Vec initial_guess = assembler.CreateConstantInitialGuess(1.0);
        // This problem seems unusally sensitive to the initial guess. Various other
        // choices failed to converge.
        Vec answer = assembler.Solve(initial_guess, true);
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double u = sqrt(x*(4-x));
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.001);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void TestHeatEquationWithNeumannOnUnitDisc( void )
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearLinearHeatEquationPde<2> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc(mesh.GetNumNodes());
        // du/dn = -0.5 on r=1
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_boundary_condition;
        p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(1), p_boundary_condition);
        
        SimpleNonlinearEllipticAssembler<2,2> assembler(&mesh, &pde, &bcc);
        
        // Set up initial Guess
        Vec initial_guess = assembler.CreateConstantInitialGuess(1.0);
        
        Vec answer = assembler.Solve(initial_guess, true);
        
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            c_vector<double, 2> r;
            r(0) = mesh.GetNodeAt(global_index)->GetPoint()[0];
            r(1) = mesh.GetNodeAt(global_index)->GetPoint()[1];
            double u = -0.25 * inner_prod(r, r) + 2.25;
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.01);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void TestWithHeatEquation2DAndNeumannBCs()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearHeatEquationPde<2> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc(mesh.GetNumNodes());
        // u(y=0) = 0
        ConstBoundaryCondition<2>* zero_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
        while (node_iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            if (fabs((*node_iter)->GetPoint()[1]) < 1e-12)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, zero_boundary_condition);
            }
            node_iter++;
        }
        
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
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
        
        Vec answer = assembler.Solve(initial_guess, true);
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
            double u = sqrt(y*(4-y));
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.15);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void Test2dOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        Example2DNonlinearEllipticPde pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<2>* p_boundary_condition;
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
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
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator elt_iter = mesh.GetBoundaryElementIteratorBegin();
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
        Vec answer = assembler.Solve(initial_guess);
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
            double u = 1 + x*x + y*y;
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.01);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(answer);
        
        // Analytical Jacobian
        answer=assembler.Solve(initial_guess, true);
        
        // Check result
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
            double u = 1 + x*x + y*y;
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.01);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
    void TestNasty2dEquationOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        ExampleNasty2dNonlinearEllipticPde pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<2>* p_boundary_condition;
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetBoundaryNodeIteratorBegin();
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
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator elt_iter = mesh.GetBoundaryElementIteratorBegin();
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
        Vec answer = assembler.Solve(initial_guess);
        
        // Check result
        double *p_answer;
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
            double u = 1 + sin(x)*sin(x) + y*y;
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.01);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(answer);
        
        // Analytical Jacobian
        answer=assembler.Solve(initial_guess, true);
        
        // Check result
        VecGetArray(answer, &p_answer);
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
            double u = 1 + sin(x)*sin(x) + y*y;
            TS_ASSERT_DELTA(p_answer[local_index], u, 0.01);
        }
        VecRestoreArray(answer, &p_answer);
        VecDestroy(initial_guess);
        VecDestroy(answer);
    }
    
};

#endif //_TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
