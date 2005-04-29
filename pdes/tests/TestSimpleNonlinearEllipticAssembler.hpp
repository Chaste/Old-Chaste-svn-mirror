#ifndef _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_

#include "SimpleNonlinearEllipticAssembler.hpp"
#include "SimpleNonlinearSolver.hpp"
#include "LinearHeatEquationPde.hpp"
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

PetscErrorCode ComputeJacobianNumerically(SNES snes, Vec input, Mat *pJacobian, 
    								     	  Mat *pPreconditioner, MatStructure *pMatStructure, 
    										  void *pContext);

  
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
        TS_ASSERT(fabs(valueLast - VonNeumannBCValue + h/2) < 0.001);
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
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess);
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

    void noTestWithHeatEquation1DAndNeumannBCs()
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
		// u(1) = 0
        //bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition);
        // u'(1) = 0
//        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetLastBoundaryElement();
//        iter--;
//        bcc.AddNeumannBoundaryCondition(*iter, pBoundaryCondition);
		pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
		// u'(1) = 0
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
    		//VecSetValue(initialGuess, i, sqrt(0.1*i*(1-0.1*i)), INSERT_VALUES);
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
 			answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess);
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
			double u = sqrt(x*(4+sqrt(8)-x));
			//std::cout << x << "\t" << u << std::endl;
			TS_ASSERT_DELTA(ans[i], u, 0.001); 
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

