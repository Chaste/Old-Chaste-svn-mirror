#ifndef _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_

#include "SimpleNonlinearEllipticAssembler.hpp"
#include <cxxtest/TestSuite.h>
//#include "petscvec.h"
//#include "petscmat.h"
  
class TestSimpleNonlinearEllipticAssembler : public CxxTest::TestSuite 
{
	
public:
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
				
        /*AbstractLinearEllipticPde<SPACE_DIM> *pPde; */
        // Instantiate PDE object
		//LinearHeatEquationPde<1> pPde;  
		
		
		// Boundary conditions
        BoundaryConditionsContainer<1,1> rBoundaryConditions;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
        rBoundaryConditions.AddDirichletBoundaryCondition(rMesh.GetNodeAt(0), pBoundaryCondition);
        
                       
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
     	 	     	
     	TS_ASSERT_THROWS_NOTHING(Vec Result = ComputeResidual(rMesh,
                       /*pPde, */
                       rBoundaryConditions,                       
                       currentSolution_vector/*,
                       pGaussianQuadratureRule*/));
       
       // Test that if we pass in a vector of ones, then in each element 
       // residual should return -h (i.e. -0.15)
        for (int i = 0; i<rMesh.GetNumNodes(); i++)
 		{
     		VecSetValue(currentSolution_vector, i, (PetscReal) 1, INSERT_VALUES);
 		}
       Vec Result = ComputeResidual(rMesh, rBoundaryConditions, currentSolution_vector);
        
                   
        PetscScalar *answerElements;
        VecGetArray(Result, &answerElements);
        double value1 = answerElements[1];
        double value2 = answerElements[2];     
		VecRestoreArray(Result,&answerElements);   
//		std::cout<< "Residual: 1st entry is "  << value1 << std::endl;       
//		std::cout<< "Residual: 2nd entry is "  << value2 << std::endl;    
        TS_ASSERT(fabs(value1 + 0.15) < 0.001);
        TS_ASSERT(fabs(value2 + 0.15) < 0.001);
        
                       
       
     }
        
        
        
};

#endif //_TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
