#ifndef _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "petscmat.h"
#include <vector>
#include <cmath>
#include <iostream>

#include "SimpleNonlinearEllipticAssembler.hpp"
#include "SimpleNonlinearSolver.hpp"

#include "Node.hpp" 
#include "Element.hpp"
#include "ConformingTetrahedralMesh.cpp"

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
#include "TrianglesMeshReader.hpp"

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
	Vec CreateConstantInitialGuessVec(int size, double value)
	{
		Vec initial_guess = CreateInitialGuessVec(size);
		VecSet(&value, initial_guess);
		VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess);
		return initial_guess;
	}

public:
		
	void TestComputeResidual( void )
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/practical1_1d_mesh"); 
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Boundary conditions
		BoundaryConditionsContainer<1,1> boundary_conditions(1, mesh.GetNumNodes());
		//Adding Dirichlet BC at node 0
		double DirichletBCValue = 5.0;
		ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(DirichletBCValue);
		boundary_conditions.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);

		// adding von Neumann BC at the last node
		double VonNeumannBCValue = 9.0;
		ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(VonNeumannBCValue);
		ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetLastBoundaryElement();
		iter--; // to be consistent with c++ :))), GetLastBoundaryElement points to one element passed it
		boundary_conditions.AddNeumannBoundaryCondition(*iter,pBoundaryCondition1);

		// initialize currentSolution_vector
		Vec currentSolution_vector;
		VecCreate(PETSC_COMM_WORLD, &currentSolution_vector);
		VecSetSizes(currentSolution_vector,PETSC_DECIDE,mesh.GetNumNodes());
		VecSetFromOptions(currentSolution_vector);
		
		Vec res_vector;
		VecDuplicate(currentSolution_vector,&res_vector);
		
		SimpleNonlinearEllipticAssembler<1,1> assembler;
		
		assembler.mpMesh = &mesh;
		NonlinearHeatEquationPde<1> pde;
		assembler.mpPde = &pde;
		assembler.mpBoundaryConditions = &boundary_conditions;
	
		assembler.ComputeResidual(currentSolution_vector, res_vector);
	
		// Set current solution to 1 and compute residual
		double h = 0.01;
		for (int i = 0; i<mesh.GetNumNodes(); i++)
		{
			VecSetValue(currentSolution_vector, i, (PetscReal) 1, INSERT_VALUES);
		}
		double InitialGuess = 1.0;
		Vec result;
		VecDuplicate(currentSolution_vector, &result);

		//15-SEP-2005 This is where we got to....
        assembler.ComputeResidual(currentSolution_vector, result);
		
 		PetscScalar *answerElements;
		VecGetArray(result, &answerElements);

        int lo, hi;
        VecGetOwnershipRange(result,&lo,&hi);
  //      VecView(result,     PETSC_VIEWER_STDOUT_WORLD );
  
        if (lo<=0 && 0<hi)
        {
		    double value1 = answerElements[0-lo];
            TS_ASSERT(fabs(value1 + DirichletBCValue - InitialGuess) < 0.001);
        }
        if (lo<=1 && 1<hi)
		{
            double value2 = answerElements[1-lo];
            TS_ASSERT(fabs(value2 + h) < 0.001);
        }
        if (lo<=mesh.GetNumNodes()-1 && mesh.GetNumNodes()-1<hi)
		{
            double valueLast = answerElements[mesh.GetNumNodes()-1-lo];
		    TS_ASSERT(fabs(valueLast + VonNeumannBCValue + h/2) < 0.001);
        }
        VecRestoreArray(result, &answerElements);   
		VecDestroy(result);
		VecDestroy(res_vector);
		VecDestroy(currentSolution_vector);
	}

	/**
	 * \todo This should be made into a proper test that the jacobian calculation
	 * is correct, for some test cases.
	 */
	void TestComputeJacobianNumerically(void)
	{	
		SNES snes;
		
		// Set up input vector - not actually used but need to be passed!
		Vec input;
		VecCreate(PETSC_COMM_WORLD, &input);
		VecSetSizes(input,PETSC_DECIDE,2);
		VecSetFromOptions(input);
		VecSetValue(input, 0, (PetscReal) 0, INSERT_VALUES);
		VecSetValue(input, 1, (PetscReal) 0, INSERT_VALUES);
		   
		// Set up Jacobian matrix - results written into this Mat object
		Mat jacobian;
		MatCreate(PETSC_COMM_WORLD, PETSC_DETERMINE, PETSC_DETERMINE, 2, 2, &jacobian);
		MatSetFromOptions(jacobian);

		//int errcode = ComputeJacobianPetsc(snes, input, &jacobian, NULL, NULL, NULL);

		//std::cout << "Our J matrix: " << std::endl;
		//MatView(jacobian,0);
		
		VecDestroy(input);
		MatDestroy(jacobian);
	}
	
	void TestNumericalAgainstAnalyticJacobian()
	{
		Mat numerical_jacobian;
		MatCreate(PETSC_COMM_WORLD, PETSC_DETERMINE, PETSC_DETERMINE, 11, 11, &numerical_jacobian);
		MatSetType(numerical_jacobian, MATSEQDENSE);
	
		Mat analytic_jacobian;
		MatCreate(PETSC_COMM_WORLD, PETSC_DETERMINE, PETSC_DETERMINE, 11, 11, &analytic_jacobian);
		MatSetType(analytic_jacobian, MATSEQDENSE);
	
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		// Instantiate PDE object
		NonlinearHeatEquationPde<1> pde;  
		
		// Boundary conditions
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
		//pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition);
			
		SimpleNonlinearEllipticAssembler<1,1> assembler;
			
		// Set up initial solution guess for residuals
		int length=mesh.GetNumNodes();
		Vec initial_guess;
		VecCreate(PETSC_COMM_WORLD, &initial_guess);
		VecSetSizes(initial_guess, PETSC_DECIDE, length);
		VecSetType(initial_guess, VECSEQ);
		for(int i=0; i<length ; i++)
		{
			//VecSetValue(initial_guess, i, sqrt(0.1*i*(1-0.1*i)), INSERT_VALUES);
			//VecSetValue(initial_guess, i, 0.25, INSERT_VALUES);
			VecSetValue(initial_guess, i, (-0.01*i*i), INSERT_VALUES);
		}
		VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess); 
			
		// Store data structures as object members
		assembler.mpMesh = &mesh;
		assembler.mpPde = &pde;
		assembler.mpBoundaryConditions = &bcc;
		
		int errcode = assembler.ComputeJacobianNumerically(initial_guess, &numerical_jacobian);
		TS_ASSERT_EQUALS(errcode, 0);
									
		errcode = assembler.ComputeJacobianAnalytically(initial_guess, &analytic_jacobian);
		TS_ASSERT_EQUALS(errcode, 0);
	
		MatAssemblyBegin(numerical_jacobian,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(numerical_jacobian,MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(analytic_jacobian,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(analytic_jacobian,MAT_FINAL_ASSEMBLY);
	
//		TS_TRACE("Numerical:");
//		MatView(numerical_jacobian, 0);
//		TS_TRACE("Analytical:");
//		MatView(analytic_jacobian, 0);
			
		PetscScalar numerical[11*11], analytic[11*11];
		PetscInt ids[11], n=11;
		for (int i=0; i<n; i++)
		{
			ids[i] = i;
		}
		
		// Check matrices are the same, to within numerical error tolerance
		MatGetValues(numerical_jacobian,n,ids,n,ids,numerical);
		MatGetValues(analytic_jacobian,n,ids,n,ids,analytic);
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
				TS_ASSERT_DELTA(numerical[i*n+j], analytic[i*n+j], 0.001);
			}
		}

		VecDestroy(initial_guess);
		MatDestroy(numerical_jacobian);
		MatDestroy(analytic_jacobian);
	}

	void TestWithHeatEquation1D()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquationPde<1> pde;  
		
		// Boundary conditions
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition);
		
		SimpleNonlinearEllipticAssembler<1,1> assembler;
		SimpleNonlinearSolver solver;
		
		// Set up initial guess
		Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
		for (int i=0; i<mesh.GetNumNodes(); i++)
		{
			VecSetValue(initial_guess, i, (-0.01*i*i), INSERT_VALUES);
		}
		VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess); 

		Vec answer;

		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
			 
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = sqrt(x*(1-x));
			TS_ASSERT_DELTA(ans[i], u, 0.001); 
		} 
		ierr = VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}

	void TestWithHeatEquation1DAndNeumannBCs()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquationPde<1> pde;
		 
		// Boundary conditions
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		// u(0) = 0
		ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
		// u(1)*u'(1) = 1
		pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
		ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetLastBoundaryElement();
		iter--;
		bcc.AddNeumannBoundaryCondition(*iter, pBoundaryCondition);

		// Nonlinear solver to use
		SimpleNonlinearSolver solver;
		
		// Set up initial Guess
		Vec initial_guess = CreateConstantInitialGuessVec(mesh.GetNumNodes(), 0.25);
		
		// Nonlinear assembler to use
		SimpleNonlinearEllipticAssembler<1,1> assembler;
		
		// Set no. of gauss points to use
		assembler.SetNumberOfQuadraturePointsPerDimension(3);

		// Solve the PDE
		Vec answer;
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
				
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = sqrt(x*(4-x));
			TS_ASSERT_DELTA(ans[i], u, 0.001);
		}
		VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}

	void TestWithHeatEquation1D2()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquation2Pde<1> pde;  
		
		// Boundary conditions
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
		ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(exp(1.0));
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition1);
		
		SimpleNonlinearEllipticAssembler<1,1> assembler;
		SimpleNonlinearSolver solver;
		
		// Set up initial Guess
		Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
		for (int i=0; i<mesh.GetNumNodes(); i++)
		{
			VecSetValue(initial_guess, i, (1.0+0.01*i*i), INSERT_VALUES);
		}
		VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess); 
				
		Vec answer;
		
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}

		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = exp(0.5*(3.0*x-x*x));
			//std::cout << x << "\t" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.001); 
		}
		
		VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}
	
	void TestWithHeatEquation1D3()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquation3Pde<1> pde;  
		
		// Boundary conditions
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(sqrt(2.0));
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
		ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0); 
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition1);
		

		SimpleNonlinearEllipticAssembler<1,1> assembler(3);
		SimpleNonlinearSolver solver;
		 
		// Set up initial Guess
		Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
		for (int i=0; i<mesh.GetNumNodes(); i++)
		{
			VecSetValue(initial_guess, i, (1.5-0.15*i), INSERT_VALUES);
		}
		VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess); 

		Vec answer;
		
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
				
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = sqrt(2.0*(exp(-x)-x*exp(-1.0)));
			TS_ASSERT_DELTA(ans[i], u, 0.001); 
		} 
		ierr = VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}
	
	void TestWithHeatEquation1D4()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquation4Pde<1> pde;  
		
		// Boundary conditions
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		// u(1) = exp(1.0)
		ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(exp(-1.0));
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition);
		// u(0)^2*u'(0) = 0.0
		pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
		ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
		bcc.AddNeumannBoundaryCondition(*iter, pBoundaryCondition);

		SimpleNonlinearEllipticAssembler<1,1> assembler;
		SimpleNonlinearSolver solver;
		 
		// Set up initial Guess
		Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
		double x1;
		for (int i=0; i<mesh.GetNumNodes(); i++)
		{
			x1=0.1*(double)(i);
			VecSetValue(initial_guess, i, 0.35*(1-x1*x1), INSERT_VALUES);
		}
		VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess); 

		Vec answer;
		
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
				
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = x*exp(-x);
			TS_ASSERT_DELTA(ans[i], u, 0.01); 
		} 
		VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}

	void TestWithHeatEquation1D5()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquation5Pde<1> pde;  

		// Boundary conditions
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		// u(1) = exp(-1.0)
		ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(exp(-1.0));
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition);
		// u(0)^2*u'(0) = -1.0
		// Note that we specify 1 as the value, since figuring out which direction
		// the normal is in is hard in 1D.
		pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
		ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
		bcc.AddNeumannBoundaryCondition(*iter, pBoundaryCondition);

		SimpleNonlinearEllipticAssembler<1,1> assembler;
		SimpleNonlinearSolver solver;
		 
		// Set up initial Guess
		Vec initial_guess = CreateInitialGuessVec(mesh.GetNumNodes());
		double x1;
		for(int i=0; i<mesh.GetNumNodes(); i++)
		{
			x1=0.1*(double)(i);
			VecSetValue(initial_guess, i, 0.35*(1-x1*x1), INSERT_VALUES);
		}
		VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess); 

		Vec answer;
		
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
				
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = exp(-x);
			TS_ASSERT_DELTA(ans[i], u, 0.01); 
		} 
		VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}

	void TestWithHeatEquation1DAndNeumannBCs2()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquationPde<1> pde;
		
		// Boundary conditions
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		// u(1) = sqrt(3)
		ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(sqrt(3));
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition);
		// u(0)*u'(0) = 2
		// Note that we specify -2 as the value, since figuring out which direction
		// the normal is in is hard in 1D.
		pBoundaryCondition = new ConstBoundaryCondition<1>(-2.0);
		ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
		bcc.AddNeumannBoundaryCondition(*iter, pBoundaryCondition);

		SimpleNonlinearEllipticAssembler<1,1> assembler;
		SimpleNonlinearSolver solver;
		
		// Set up initial Guess
		Vec initial_guess = CreateConstantInitialGuessVec(mesh.GetNumNodes(), 1.0);
		// This problem seems unusally sensitive to the initial guess. Various other
		// choices failed to converge.
		
		Vec answer;
		
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
				
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = sqrt(x*(4-x));
			TS_ASSERT_DELTA(ans[i], u, 0.001);
		}
		VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}
	
	void TestHeatEquationWithNeumannOnUnitDisc( void )
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/disk_522_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearLinearHeatEquationPde<2> pde;
		
		// Boundary conditions
		BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
		// du/dn = -0.5 on r=1
		ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
		ConstBoundaryCondition<2>* pBoundaryCondition;
		pBoundaryCondition = new ConstBoundaryCondition<2>(-0.5);
		while (iter != mesh.GetLastBoundaryElement())
		{
			bcc.AddNeumannBoundaryCondition(*iter, pBoundaryCondition);
			iter++;
		}
		// u = 2 at some point on the boundary, say node 1
		pBoundaryCondition = new ConstBoundaryCondition<2>(2.0);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(1), pBoundaryCondition);
		
		SimpleNonlinearEllipticAssembler<2,2> assembler;
		SimpleNonlinearSolver solver;
		
		// Set up initial Guess
		Vec initial_guess = CreateConstantInitialGuessVec(mesh.GetNumNodes(), 1.0);
		
		Vec answer;
		
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
		
		double *res;
		int ierr = VecGetArray(answer, &res);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			VectorDouble r(2);
			r(0) = mesh.GetNodeAt(i)->GetPoint()[0];
			r(1) = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = -0.25 * r.L2Norm() * r.L2Norm() + 2.25;
			TS_ASSERT_DELTA(res[i], u, 0.01);
		}
		VecRestoreArray(answer, &res);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}
	
	void TestWithHeatEquation2DAndNeumannBCs()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquationPde<2> pde;
		 
		// Boundary conditions
		BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
		// u(y=0) = 0
		ConstBoundaryCondition<2>* zeroBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
		ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetFirstBoundaryNode();
		while (node_iter != mesh.GetLastBoundaryNode())
		{
			if (fabs((*node_iter)->GetPoint()[1]) < 1e-12)
			{
				bcc.AddDirichletBoundaryCondition(*node_iter, zeroBoundaryCondition);
			}
			node_iter++;
		}

		ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
		FunctionalBoundaryCondition<2>* oneBoundaryCondition = new FunctionalBoundaryCondition<2>(&one_bc);
		AbstractBoundaryCondition<2>* pBoundaryCondition;
		while (iter != mesh.GetLastBoundaryElement())
		{
			double x = (*iter)->GetNodeLocation(0,0);
			double y = (*iter)->GetNodeLocation(0,1);
			if (fabs(y-1.0) < 1e-12)
			{
				// u(y=1)*u'(y=1) = 1
				pBoundaryCondition = oneBoundaryCondition;
			}
			else
			{
				// No flux across left & right
				pBoundaryCondition = zeroBoundaryCondition;
			}
				
			bcc.AddNeumannBoundaryCondition(*iter, pBoundaryCondition);
				
			iter++;
		}
	
		SimpleNonlinearEllipticAssembler<2,2> assembler;
		SimpleNonlinearSolver solver;
			
		// Set up initial Guess
   		Vec initial_guess = CreateConstantInitialGuessVec(mesh.GetNumNodes(), 0.25);
		
		Vec answer;
		
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
				
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = sqrt(y*(4-y));
			TS_ASSERT_DELTA(ans[i], u, 0.15);
		}
		VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}

	void Test2dOnUnitSquare()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);

		// Instantiate PDE object
		Example2DNonlinearEllipticPde pde;

		// Boundary conditions
		BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
		ConstBoundaryCondition<2>* pBoundaryCondition;
		ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetFirstBoundaryNode();
		while (node_iter != mesh.GetLastBoundaryNode())
		{
			double x = (*node_iter)->GetPoint()[0];
			double y = (*node_iter)->GetPoint()[1];
			pBoundaryCondition = NULL;
			if (fabs(x) < 1e-12)
			{
				// On x=0, u=1+y^2
				pBoundaryCondition = new ConstBoundaryCondition<2>(1 + y*y);
			}
			else if (fabs(y) < 1e-12)
			{
				// On y=0, u=1+x^2
				pBoundaryCondition = new ConstBoundaryCondition<2>(1 + x*x);
			}
			if (pBoundaryCondition)
			{
				bcc.AddDirichletBoundaryCondition(*node_iter, pBoundaryCondition);
			}
				
			node_iter++;
		}
		FunctionalBoundaryCondition<2>* pBC;
		ConformingTetrahedralMesh<2,2>::BoundaryElementIterator elt_iter = mesh.GetFirstBoundaryElement();
		while (elt_iter != mesh.GetLastBoundaryElement())
		{
			double x = (*elt_iter)->GetNodeLocation(0,0);
			double y = (*elt_iter)->GetNodeLocation(0,1);
			pBC = NULL;
			if (fabs(y-1.0) < 1e-12)
			{
				// On y=1, Dgradu_dot_n = 2(2+x^2)
				pBC = new FunctionalBoundaryCondition<2>(&bc_y1_func);
			}
			else if (fabs(x-1.0) < 1e-12)
			{
				// On x=1, Dgradu_dot_n = 2(2+y^2)
				pBC = new FunctionalBoundaryCondition<2>(&bc_x1_func);
			}
			if (pBC)
			{
				bcc.AddNeumannBoundaryCondition(*elt_iter, pBC);
			}
			
			elt_iter++;
		}

		SimpleNonlinearEllipticAssembler<2,2> assembler;
		SimpleNonlinearSolver solver;
		
		// Set up initial Guess
		Vec initial_guess = CreateConstantInitialGuessVec(mesh.GetNumNodes(), 4.0);
		
		Vec answer;
		
		// Numerical Jacobian
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
		
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = 1 + x*x + y*y;
			TS_ASSERT_DELTA(ans[i], u, 0.01);
		}
		VecRestoreArray(answer, &ans);
		VecDestroy(answer);
		
		// Analytical Jacobian
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
		
		// Check result
		ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = 1 + x*x + y*y;
			TS_ASSERT_DELTA(ans[i], u, 0.01);
		}
		VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}

	void TestNasty2dEquationOnUnitSquare()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);

		// Instantiate PDE object
		ExampleNasty2dNonlinearEllipticPde pde; 

		// Boundary conditions
		BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
		ConstBoundaryCondition<2>* pBoundaryCondition;
		ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetFirstBoundaryNode();
		while (node_iter != mesh.GetLastBoundaryNode())
		{
			double x = (*node_iter)->GetPoint()[0];
			double y = (*node_iter)->GetPoint()[1];
			pBoundaryCondition = NULL;
			if (fabs(x) < 1e-12)
			{
				// On x=0, u=1+y^2
				pBoundaryCondition = new ConstBoundaryCondition<2>(1 + y*y);
			}
			else if (fabs(y) < 1e-12)
			{
				// On y=0, u=1+sin^2(x)
				pBoundaryCondition = new ConstBoundaryCondition<2>(1 + sin(x)*sin(x));
			}
			if (pBoundaryCondition)
			{
				bcc.AddDirichletBoundaryCondition(*node_iter, pBoundaryCondition);
			}
			
			node_iter++;
		}
		FunctionalBoundaryCondition<2>* pBC;
		ConformingTetrahedralMesh<2,2>::BoundaryElementIterator elt_iter = mesh.GetFirstBoundaryElement();
		while (elt_iter != mesh.GetLastBoundaryElement())
		{
			double x = (*elt_iter)->GetNodeLocation(0,0);
			double y = (*elt_iter)->GetNodeLocation(0,1);
			pBC = NULL;
			if (fabs(y-1.0) < 1e-12)
			{
				// On y=1, Dgradu_dot_n = 2(2+sin^2(x))
				pBC = new FunctionalBoundaryCondition<2>(&bc_y1_func2);
			}
			else if (fabs(x-1.0) < 1e-12)
			{
				// On x=1, Dgradu_dot_n = sin(2)(sin^2(1)+1+y^2)
				pBC = new FunctionalBoundaryCondition<2>(&bc_x1_func2);
			}
			if (pBC)
			{
				bcc.AddNeumannBoundaryCondition(*elt_iter, pBC);
			}
			
			elt_iter++;
		}

		SimpleNonlinearEllipticAssembler<2,2> assembler;
		SimpleNonlinearSolver solver;
		
		// Set up initial Guess
		Vec initial_guess = CreateConstantInitialGuessVec(mesh.GetNumNodes(), 4.0);
		
		Vec answer;
		
		// Numerical Jacobian
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
		
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = 1 + sin(x)*sin(x) + y*y;
			TS_ASSERT_DELTA(ans[i], u, 0.01);
		}
		VecRestoreArray(answer, &ans);
		VecDestroy(answer);
		
		// Analytical Jacobian
		try {
			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initial_guess, true);
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
			TS_ASSERT(0);
		}
		
		// Check result
		ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = 1 + sin(x)*sin(x) + y*y;
			TS_ASSERT_DELTA(ans[i], u, 0.01);
		}
		VecRestoreArray(answer, &ans);
		VecDestroy(initial_guess);
		VecDestroy(answer);
	}

};

#endif //_TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
