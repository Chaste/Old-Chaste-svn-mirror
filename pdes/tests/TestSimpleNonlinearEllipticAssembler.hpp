#ifndef _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_

#include "SimpleNonlinearEllipticAssembler.hpp"
#include "SimpleNonlinearSolver.hpp"

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "petscmat.h"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <cmath>
#include <iostream>
#include "Node.hpp" 
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "FunctionalBoundaryCondition.hpp"

#include "NonlinearHeatEquationPde.hpp"
#include "NonlinearHeatEquation2Pde.hpp"
#include "NonlinearHeatEquation3Pde.hpp"
#include "Example2DNonlinearEllipticPde.hpp"
#include "NonlinearLinearHeatEquationPde.hpp"
#include "ExampleNasty2dNonlinearEllipticPde.hpp"


PetscErrorCode ComputeJacobianNumerically(SNES snes, Vec input, Mat *pJacobian, 
    								     	  Mat *pPreconditioner, MatStructure *pMatStructure, 
    										  void *pContext);

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
	
public:
	void setUp()
    {
		int FakeArgc=0;
		char *FakeArgv0="testrunner";
		char **FakeArgv=&FakeArgv0;
	    	
		PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
    }

    void testASimpleNonlinearEllipticAssembler( void )
    {
        //create a new SimpleNonlinearEllipticAssembler
        SimpleNonlinearEllipticAssembler<1,1> assembler;
    }
        
    void DONOTtestSimpleNonlinearEllipticAssembler( void )
    {
        //create a new SimpleNonlinearEllipticAssembler
        SimpleNonlinearEllipticAssembler<1,1> assembler;
//     
//        Element<1,1> element(nodes);
//        
//        AbstractBasisFunction<SPACE_DIM> *pBasisFunction;
//        pBasisFunction = new LinearBasisFunction<1>():
//
//		
//		
//		ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
//                       AbstractNonlinearEllipticPde<SPACE_DIM> *pPde, 
//                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
//                       AbstractNonlinearSolver *pSolver,
//                       GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule);
//                       
//                       
//		assembler = AssembleSystem(rMesh, pPde, rBoundaryConditions, pSolver, pBasisFunction, pGaussianQuadratureRule);
    }
        
        
     void testComputeResidual( void )
     {
		// Create mesh from mesh reader
		//TrianglesMeshReader mesh_reader("pdes/tests/meshdata/trivial_1d_mesh"); 
        //double h = 0.15;
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh"); 
        double h = 0.01;
        
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
				
		// Boundary conditions
        BoundaryConditionsContainer<1,1> boundary_conditions;
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
     	VecSetType(currentSolution_vector, VECSEQ);
     	
     	Vec res_vector;
    	VecDuplicate(currentSolution_vector,&res_vector);
        
        /* GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule;*/
     	
     	SNES snes;
     	SimpleNonlinearEllipticAssembler<1,1> assembler;
     		     	
     	assembler.mpMesh = &mesh;
     	NonlinearHeatEquationPde<1> pde;
		assembler.mpPde = &pde;
		assembler.mpBoundaryConditions = &boundary_conditions;
		LinearBasisFunction<1> basis_function;
		assembler.mpBasisFunction = &basis_function;
		//assembler.mpGaussianQuadratureRule=pGaussianQuadratureRule;

     	//TS_ASSERT_THROWS_NOTHING(
     	ComputeResidual<1, 1>(snes, currentSolution_vector, res_vector, &assembler);


		// Set current solution to 1 and compute residual
        for (int i = 0; i<mesh.GetNumNodes(); i++)
 		{
     		VecSetValue(currentSolution_vector, i, (PetscReal) 1, INSERT_VALUES);
 		}
		double InitialGuess = 1.0;
 		Vec Result;
 		VecDuplicate(currentSolution_vector, &Result);
        ComputeResidual<1,1>(snes, currentSolution_vector, Result, &assembler);
        
        PetscScalar *answerElements;
        VecGetArray(Result, &answerElements);
        double value1 = answerElements[0];
        double value2 = answerElements[1];
        double valueLast = answerElements[mesh.GetNumNodes()-1];
		VecRestoreArray(Result,&answerElements);   

        TS_ASSERT(fabs(value1 + DirichletBCValue - InitialGuess) < 0.001);
        TS_ASSERT(fabs(value2 + h) < 0.001);
        TS_ASSERT(fabs(valueLast + VonNeumannBCValue + h/2) < 0.001);
	}

     /**
     * Function ComputeJacobianNumerically() is modified and appended towards bottom
     * of this file. Function modified and doesn't actually call ComputeResidual()
     * but 'residual' and 'perturbedResidual' are hardcoded.
     * 'num_nodes' hardcoded to 2.
     */
    void testComputeJacobianNumerically(void)
    {    
    	SNES snes;
    	
    	// Set up input vector - not actually used but need to be passed!
    	Vec input;
    	VecCreate(PETSC_COMM_WORLD, &input);
		VecSetSizes(input,PETSC_DECIDE,2);
		//VecSetType(input, VECSEQ);
	   	VecSetFromOptions(input);
	   	VecSetValue(input, 0, (PetscReal) 0, INSERT_VALUES);
	   	VecSetValue(input, 1, (PetscReal) 0, INSERT_VALUES);
	   	
	   	// Set up Jacobian matrix - results written into this Mat object
   		Mat jacobian;
   		MatCreate(PETSC_COMM_WORLD, 2, 2, PETSC_DETERMINE, PETSC_DETERMINE, &jacobian);
   		MatSetType(jacobian, MATSEQDENSE);
   		//MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2, 2, &jacobian);
    	//MatSetType(jacobian, MATMPIDENSE);

    	
        int errcode = ComputeJacobianNumerically(snes, input, &jacobian, NULL, NULL, NULL);

        //std::cout << "Our J matrix: " << std::endl;
        //MatView(pJacobian,0);
        
        TS_ASSERT(true); 
    }   
    
    
    void testWithHeatEquation1D()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquationPde<1> pde;  
		
		// Boundary conditions
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        //pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition);
        
		SimpleNonlinearEllipticAssembler<1,1> assembler;
    	SimpleNonlinearSolver solver;
    	
    	// Set up solution guess for residuals
    	int length=mesh.GetNumNodes();
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		//VecSetValue(initialGuess, i, sqrt(0.1*i*(1-0.1*i)), INSERT_VALUES);
    		//VecSetValue(initialGuess, i, 0.25, INSERT_VALUES);
    		VecSetValue(initialGuess, i, (-0.01*i*i), INSERT_VALUES);
    	}
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess); 
		
		//
		GaussianQuadratureRule<1> quadRule(2);
		LinearBasisFunction<1> basis_func;
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	
    	//TS_TRACE("Calling AssembleSystem");
 		try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
    	//TS_TRACE("System solved");
    	 	
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = sqrt(x*(1-x));
			//std::cout << x << "\t" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.001); 
		} 
		VecRestoreArray(answer, &ans);
	}
	
	void TestNumericalAgainstAnalyticJacobian()
	{
		TestStuff();
	}

    void TestWithHeatEquation1DAndNeumannBCs()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquationPde<1> pde;
		 
		// Boundary conditions
        BoundaryConditionsContainer<1,1> bcc;
        // u(0) = 0
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
		// u(1)*u'(1) = 1
		pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetLastBoundaryElement();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, pBoundaryCondition);

		SimpleNonlinearEllipticAssembler<1,1> assembler;
    	SimpleNonlinearSolver solver;
    	
    	// Set up solution guess for residuals
    	int length=mesh.GetNumNodes();
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		//VecSetValue(initialGuess, i, sqrt(0.1*i*(4-0.1*i)), INSERT_VALUES);
    		VecSetValue(initialGuess, i, 0.25, INSERT_VALUES);
    		//VecSetValue(initialGuess, i, (-0.01*i*i), INSERT_VALUES);
    	}
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess); 
		
		GaussianQuadratureRule<1> quadRule(2);
		LinearBasisFunction<1> basis_func;
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	
    	//TS_TRACE("Calling AssembleSystem");
    	try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
 		//TS_TRACE("System solved");
    	    	
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = sqrt(x*(4-x));
			//std::cout << x << "\t" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.001);
		}
		VecRestoreArray(answer, &ans);
	}

	void testWithHeatEquation1D2()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquation2Pde<1> pde;  
		
		// Boundary conditions
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(1.0); 
        pBoundaryCondition1 = new ConstBoundaryCondition<1>(exp(1.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition1);
        
		SimpleNonlinearEllipticAssembler<1,1> assembler;
    	SimpleNonlinearSolver solver;
    	
    	// Set up solution guess for residuals
    	int length=mesh.GetNumNodes();
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		//VecSetValue(initialGuess, i, sqrt(0.1*i*(1-0.1*i)), INSERT_VALUES);
    		//VecSetValue(initialGuess, i, 0.25, INSERT_VALUES);
    		VecSetValue(initialGuess, i, (1.0+0.01*i*i), INSERT_VALUES);
    	}
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess); 
		
		//
		GaussianQuadratureRule<1> quadRule(2);
		LinearBasisFunction<1> basis_func;
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	
    	//TS_TRACE("Calling AssembleSystem");
 		try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
    	//TS_TRACE("System solved");
    	    	
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
	}
	
	void testWithHeatEquation1D3()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquation3Pde<1> pde;  
		
		// Boundary conditions
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(sqrt(2.0));
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0); 
        pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition1);
        

		SimpleNonlinearEllipticAssembler<1,1> assembler;
    	SimpleNonlinearSolver solver;
    	 
    	// Set up solution guess for residuals
    	int length=mesh.GetNumNodes();
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		//VecSetValue(initialGuess, i, sqrt(0.1*i*(1-0.1*i)), INSERT_VALUES);
    		//VecSetValue(initialGuess, i, 0.25, INSERT_VALUES);
    		VecSetValue(initialGuess, i, (1.5-0.15*i), INSERT_VALUES);
    	}
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess); 
		
		//
		GaussianQuadratureRule<1> quadRule(2);
		LinearBasisFunction<1> basis_func;
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	
    	//TS_TRACE("Calling AssembleSystem");
 		try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
    	//TS_TRACE("System solved");
    	    	
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = sqrt(2.0*(exp(-x)-x*exp(-1.0)));
			//std::cout << x << "\t" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.001); 
		} 
		VecRestoreArray(answer, &ans);
	}
	void TestWithHeatEquation1DAndNeumannBCs2()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquationPde<1> pde;
		
		// Boundary conditions
        BoundaryConditionsContainer<1,1> bcc;
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
    	
    	// Set up solution guess for residuals
    	int length=mesh.GetNumNodes();
		
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		// This problem seems unusally sensitive to the initial guess.
    		// The commented out choices (except for the right answer) failed to converge
    		//VecSetValue(initialGuess, i, sqrt(0.1*i*(4-0.1*i)), INSERT_VALUES);
    		VecSetValue(initialGuess, i, 1.0, INSERT_VALUES);
    		//VecSetValue(initialGuess, i, 0.25, INSERT_VALUES);
    		//VecSetValue(initialGuess, i, (+0.01*i*i), INSERT_VALUES);
    		//VecSetValue(initialGuess, i, 0.1*i, INSERT_VALUES);
    	}
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess); 
		
		GaussianQuadratureRule<1> quadRule(2);
		LinearBasisFunction<1> basis_func;

    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	
    	//TS_TRACE("Calling AssembleSystem");
    	try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
 		//TS_TRACE("System solved");
    	    	
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = sqrt(x*(4-x));
			//std::cout << x << "\t" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.001);
		}
		VecRestoreArray(answer, &ans);
	}
	
	void TestHeatEquationWithNeumannOnUnitDisc( void )
    {
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/disk_522_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        NonlinearLinearHeatEquationPde<2> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2> bcc;
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
    	int length=mesh.GetNumNodes();
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		VecSetValue(initialGuess, i, 1.0, INSERT_VALUES);
    		//VecSetValue(initialGuess, i, 0.25, INSERT_VALUES);
    		//VecSetValue(initialGuess, i, (+0.01*i*i), INSERT_VALUES);
    		//VecSetValue(initialGuess, i, 0.1*i, INSERT_VALUES);
    	}
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess); 
		
		GaussianQuadratureRule<2> quadRule(2);
		LinearBasisFunction<2> basis_func;

    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	
    	//TS_TRACE("Calling AssembleSystem");
    	try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
 		//TS_TRACE("System solved");
        
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
    }
    
    void TestWithHeatEquation2DAndNeumannBCs()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/square_128_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		NonlinearHeatEquationPde<2> pde;
		 
		// Boundary conditions
        BoundaryConditionsContainer<2,2> bcc;
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
    	
    	// Set up solution guess for residuals
    	int length=mesh.GetNumNodes();
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		//VecSetValue(initialGuess, i, sqrt(0.1*i*(4-0.1*i)), INSERT_VALUES);
    		VecSetValue(initialGuess, i, 0.25, INSERT_VALUES);
    		//VecSetValue(initialGuess, i, (-0.01*i*i), INSERT_VALUES);
    	}
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess); 
		
		GaussianQuadratureRule<2> quadRule(2);
		LinearBasisFunction<2> basis_func;
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	
    	//TS_TRACE("Calling AssembleSystem");
    	try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
 		//TS_TRACE("System solved");
    	    	
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = sqrt(y*(4-y));
			//std::cout << x << "\t" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.15);
		}
		VecRestoreArray(answer, &ans);
	}

    void Test2dOnUnitSquare()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/square_128_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);

		// Instantiate PDE object
		Example2DNonlinearEllipticPde pde;

		// Boundary conditions
        BoundaryConditionsContainer<2,2> bcc;
        ConstBoundaryCondition<2>* pBoundaryCondition;
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetFirstBoundaryNode();
        while (node_iter != mesh.GetLastBoundaryNode())
        {
        	double x = (*node_iter)->GetPoint()[0];
        	double y = (*node_iter)->GetPoint()[1];
        	pBoundaryCondition = NULL;
	        // On x=0, u=1+y^2
	        if (fabs(x) < 1e-12)
	        {
        		pBoundaryCondition = new ConstBoundaryCondition<2>(1 + y*y);
	        }
        	// On y=0, u=1+x^2
        	if (fabs(y) < 1e-12)
	        {
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
        	// On x=1, Dgradu_dot_n = 2(2+y^2)
	        if (fabs(x-1.0) < 1e-12)
	        {
        		pBC = new FunctionalBoundaryCondition<2>(&bc_x1_func);
	        }
	        // On y=1, Dgradu_dot_n = 2(2+x^2)
        	if (fabs(y-1.0) < 1e-12)
	        {
        		pBC = new FunctionalBoundaryCondition<2>(&bc_y1_func);
	        }
	        if (pBC)
	        {
        		bcc.AddNeumannBoundaryCondition(*elt_iter, pBC);
	        }
        
        	elt_iter++;
		}

		SimpleNonlinearEllipticAssembler<2,2> assembler;
    	SimpleNonlinearSolver solver;
    	
    	// Set up solution guess for residuals
    	int length=mesh.GetNumNodes();
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		VecSetValue(initialGuess, i, 4.0, INSERT_VALUES);
    	}
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess); 
		
		GaussianQuadratureRule<2> quadRule(2);
		LinearBasisFunction<2> basis_func;
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	
    	// Numerical Jacobian
    	try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
    	
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = 1 + x*x + y*y;
			//std::cout << "u(" << x << "," << y << ")=" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.01);
		}
		VecRestoreArray(answer, &ans);
		
		// Analytical Jacobian
    	try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
    	
		// Check result
		ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = 1 + x*x + y*y;
			//std::cout << "u(" << x << "," << y << ")=" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.01);
		}
		VecRestoreArray(answer, &ans);
	}

	
	void TestNasty2dEquationOnUnitSquare()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/square_128_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);

		// Instantiate PDE object
		ExampleNasty2dNonlinearEllipticPde pde;

		// Boundary conditions
        BoundaryConditionsContainer<2,2> bcc;
        ConstBoundaryCondition<2>* pBoundaryCondition;
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = mesh.GetFirstBoundaryNode();
        while (node_iter != mesh.GetLastBoundaryNode())
        {
        	double x = (*node_iter)->GetPoint()[0];
        	double y = (*node_iter)->GetPoint()[1];
        	pBoundaryCondition = NULL;
	        // On x=0, u=1+y^2
	        if (fabs(x) < 1e-12)
	        {
        		pBoundaryCondition = new ConstBoundaryCondition<2>(1 + y*y);
	        }
        	// On y=0, u=1+sin^2(x)
        	if (fabs(y) < 1e-12)
	        {
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
        	// On x=1, Dgradu_dot_n = sin(2)(sin^2(1)+1+y^2)
	        if (fabs(x-1.0) < 1e-12)
	        {
        		pBC = new FunctionalBoundaryCondition<2>(&bc_x1_func2);
	        }
	        // On y=1, Dgradu_dot_n = 2(2+sin^2(x))
        	if (fabs(y-1.0) < 1e-12)
	        {
        		pBC = new FunctionalBoundaryCondition<2>(&bc_y1_func2);
	        }
	        if (pBC)
	        {
        		bcc.AddNeumannBoundaryCondition(*elt_iter, pBC);
	        }
        
        	elt_iter++;
		}

		SimpleNonlinearEllipticAssembler<2,2> assembler;
    	SimpleNonlinearSolver solver;
    	
    	// Set up solution guess for residuals
    	int length=mesh.GetNumNodes();
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		VecSetValue(initialGuess, i, 4.0, INSERT_VALUES);
    	}
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess); 
		
		GaussianQuadratureRule<2> quadRule(2);
		LinearBasisFunction<2> basis_func;
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	
    	// Numerical Jacobian
    	try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
    	
		// Check result
		double *ans;
		int ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = 1 + sin(x)*sin(x) + y*y;
			//std::cout << "u(" << x << "," << y << ")=" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.01);
		}
		VecRestoreArray(answer, &ans);
		
		// Analytical Jacobian
    	try {
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
 		} catch (Exception e) {
 			TS_TRACE(e.getMessage());
 		}
    	
		// Check result
		ierr = VecGetArray(answer, &ans);
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double y = mesh.GetNodeAt(i)->GetPoint()[1];
			double u = 1 + sin(x)*sin(x) + y*y;
			//std::cout << "u(" << x << "," << y << ")=" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.01);
		}
		VecRestoreArray(answer, &ans);
	}
	
	
};


/**
 * Modified version of ComputeJacobianNumerically() for testing.
 * 
 * 'residual' and 'perturbedResidual' are hardcoded (do not actually call ComputeResidual()).
 * 'num_nodes' hardcoded to 2.
 */
//template <int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeJacobianNumerically(SNES snes, Vec input, Mat *pJacobian, 
    								     	  Mat *pPreconditioner, MatStructure *pMatStructure, 
    										  void *pContext)
{
	int ierr;
    Vec residual;
    Vec perturbedResidual;
    Vec result;
    
    // Commented out for testing!!
    //SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *integrator =
	//    ((SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*)pContext);    
    
    //int num_nodes = integrator->mpMesh->GetNumNodes();
    int num_nodes = 2;

    VecCreate(PETSC_COMM_WORLD, &residual);    
    VecCreate(PETSC_COMM_WORLD, &result);    
    VecCreate(PETSC_COMM_WORLD, &perturbedResidual);    
    
    VecSetSizes(residual,PETSC_DECIDE,num_nodes);
    VecSetSizes(result,PETSC_DECIDE,num_nodes);
    VecSetSizes(perturbedResidual,PETSC_DECIDE,num_nodes);
    
    //VecSetType(residual, VECSEQ);
    //VecSetType(result, VECSEQ);
    //VecSetType(perturbedResidual, VECSEQ);
    VecSetFromOptions(residual);
    VecSetFromOptions(result);
    VecSetFromOptions(perturbedResidual);

    Vec inputcopy;

    ierr = VecDuplicate(input,&inputcopy); CHKERRQ(ierr);
    ierr = VecCopy(input, inputcopy);CHKERRQ(ierr);
    
    // Hard coding residual and perturbedResidual to test since ComputeResidual() function
    // not complete!
    
	//ComputeResidual<ELEMENT_DIM, SPACE_DIM>(snes,input,residual,pContext);

	//***************************************************
    for (int row=0;row<num_nodes;row++)
    {
    	PetscScalar value = 1;
    	VecSetValue(residual, row, value, INSERT_VALUES);
    }
    //***************************************************
    
    double h = 0.00001;    
    PetscScalar subtract = -1;
    PetscScalar oneOverH = 1/h;
    
    
    for(int j = 0; j < num_nodes; j++)
    {
		PetscScalar *resultElements;
        ierr = VecSetValue(inputcopy,j,h, ADD_VALUES);CHKERRQ(ierr);
        
        //ComputeResidual<ELEMENT_DIM, SPACE_DIM>(snes, inputcopy, perturbedResidual,pContext);
        
        //*************************************************************
        for (int row=0;row<num_nodes;row++)
	   	{
    		int temp = 1;
    		if (row==j) temp += 1;
    		PetscScalar value2 = temp;
    		VecSetValue(perturbedResidual, row, value2, INSERT_VALUES);
    	}
        //*************************************************************
                
        ierr = VecWAXPY(&subtract,residual,perturbedResidual,result);CHKERRQ(ierr);
        ierr = VecScale(&oneOverH, result);CHKERRQ(ierr);
        
        ierr = VecGetArray(result,&resultElements);CHKERRQ(ierr);

        for (int i=0; i < num_nodes; i++)
        {
            ierr = MatSetValue(*pJacobian,i,j,resultElements[i],INSERT_VALUES);CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(result,&resultElements); CHKERRQ(ierr);
        
        ierr = VecSetValue(inputcopy,j,-h, ADD_VALUES); CHKERRQ(ierr);
    }
    
    VecDestroy(residual);
    VecDestroy(perturbedResidual);
    VecDestroy(result);
    VecDestroy(inputcopy);
 
    MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);
    return 0;
}

	
#endif //_TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_

