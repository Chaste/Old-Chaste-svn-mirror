#ifndef _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_

#include "SimpleNonlinearEllipticAssembler.hpp"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "petscmat.h"

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
     	
     	      
        // Create mesh (by hand!) (from TestSimpleLinearEllipticAssembler.hpp)
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/trivial_1d_mesh");
		ConformingTetrahedralMesh<1,1> rMesh;
		rMesh.ConstructFromMeshReader(mesh_reader);
		std::vector<Node<1>*> nodes;
		
		// Define PDE 
       // AbstractLinearEllipticPde<1> *pPde; 
       //	LinearHeatEquationPde<1> *pPde;  
		
		
		// Boundary conditions
		double DirichletBCValue = 5.0;
        BoundaryConditionsContainer<1,1> rBoundaryConditions;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(DirichletBCValue);
        rBoundaryConditions.AddDirichletBoundaryCondition(rMesh.GetNodeAt(0), pBoundaryCondition);
        // need to test non-zero Dirichlet BC and then vonNeumann BCs
                       
        // initialize currentSolution_vector
        Vec currentSolution_vector;
     	VecCreate(PETSC_COMM_WORLD, &currentSolution_vector);
     	VecSetSizes(currentSolution_vector,PETSC_DECIDE,rMesh.GetNumNodes());
     	VecSetType(currentSolution_vector, VECSEQ);
     	
// 		for (int i = 0; i<rMesh.GetNumNodes(); i++)
// 		{
//     		VecSetValue(currentSolution_vector, i, (PetscReal) 1, INSERT_VALUES);
// 		}
     	
     	
        
        /* GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule;*/
     	
     	
//     	// Vec ComputeResidual(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
//                       /*AbstractLinearEllipticPde<SPACE_DIM> *pPde, */
//                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
//                       Vec CurrentSolution/*,
//                       GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule*/)
     	 	     	
     	TS_ASSERT_THROWS_NOTHING(Vec Result1 = ComputeResidual(rMesh,
                       /*pPde,*/ 
                       rBoundaryConditions,                       
                       currentSolution_vector/*,
                       pGaussianQuadratureRule*/));
       
       // Test that if we pass in a vector of ones, then in each element 
       // residual should return -h (i.e. -0.15)
        for (int i = 0; i<rMesh.GetNumNodes(); i++)
 		{
     		VecSetValue(currentSolution_vector, i, (PetscReal) 1, INSERT_VALUES);
 		}
 	   double InitialGuess = 1.0;
       Vec Result = ComputeResidual(rMesh, /*pPde,*/ rBoundaryConditions, currentSolution_vector);
        
                   
        PetscScalar *answerElements;
        VecGetArray(Result, &answerElements);
        double value1 = answerElements[0];
        double value2 = answerElements[1];     
		VecRestoreArray(Result,&answerElements);   
//		std::cout<< "Residual: 1st entry is "  << value1 << std::endl;       
//		std::cout<< "Residual: 2nd entry is "  << value2 << std::endl;    
        TS_ASSERT(fabs(value1 + DirichletBCValue - InitialGuess) < 0.001);
        TS_ASSERT(fabs(value2 + 0.15) < 0.001);
        
                       
       
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
		VecSetType(input, VECSEQ);
	   	VecSetValue(input, 0, (PetscReal) 0, INSERT_VALUES);
	   	VecSetValue(input, 1, (PetscReal) 0, INSERT_VALUES);
	   	
	   	// Set up Jacobian matrix - results written into this Mat object
   		Mat pJacobian;
   		MatCreate(PETSC_COMM_WORLD, 2, 2, PETSC_DETERMINE, PETSC_DETERMINE, &pJacobian);
   		MatSetType(pJacobian, MATSEQDENSE);
    	
        int errcode = ComputeJacobianNumerically(snes, input, &pJacobian, NULL, NULL, NULL);
        
        //std::cout << "Our J matrix: " << std::endl;
        //MatView(pJacobian,0);
        
        TS_ASSERT(true);    
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
    
    VecSetType(residual, VECSEQ);
    VecSetType(result, VECSEQ);
    VecSetType(perturbedResidual, VECSEQ);
    
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
    TS_ASSERT( 1 < 2 ) ;
}

#endif //_TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_

