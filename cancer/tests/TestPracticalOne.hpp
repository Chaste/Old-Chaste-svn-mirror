#ifndef _TESTPRACTICALONE_HPP_
#define _TESTPRACTICALONE_HPP_

#include "SimpleLinearSolver.hpp"
#include "petscmat.h"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "SimpleNonlinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <iostream>
#include "Node.hpp" 
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ColumnDataWriter.hpp"

#include "Practical1Question1Pde.hpp"
#include "Practical1Question1PdeNonlinear.hpp"
#include "Practical1Question2Pde.hpp"
#include "Practical1Question3Pde.hpp"
#include "BoundaryOdeSystem.hpp"

#include "SimpleLinearEllipticAssembler.hpp"
#include "SimpleNonlinearEllipticAssembler.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "Node.hpp"

  
class TestPracticalOne : public CxxTest::TestSuite 
{
public:
	
	void testQuestion1(void)
	{
		int FakeArgc=0;
		char *FakeArgv0="testrunner";
		char **FakeArgv=&FakeArgv0;
		PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
		
		// Create mesh from mesh reader
		
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		// Instantiate PDE object
		Practical1Question1Pde<1> pde;  
		
		// Boundary conditions
		// u'(0)=0, u(1)=1
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
        bcc.AddNeumannBoundaryCondition(*iter,pBoundaryCondition1);
        
        ConstBoundaryCondition<1>* pBoundaryCondition2 = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pBoundaryCondition2);
        
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
		
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution should be u = 0.5*x*(3-x)
		for (int i=0; i < mesh.GetNumElements()+1; i++)
		{
			double x = 0.0 + 0.01*i;
			double u = 0.5*(x*x+1);
			TS_ASSERT_DELTA(res[i], u, 0.001);
		}
		VecRestoreArray(result, &res);
	}
    
    void testQuestion1Nonlinear(void)
    {
        int FakeArgc=0;
        char *FakeArgv0="testrunner";
        char **FakeArgv=&FakeArgv0;
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
        
        // Create mesh from mesh reader
        
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        // Instantiate PDE object
        Practical1Question1PdeNonlinear<1> pde;  
        
        // Boundary conditions
        // u'(0)=0, u(1)=1
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
        bcc.AddNeumannBoundaryCondition(*iter,pBoundaryCondition1);
        
        ConstBoundaryCondition<1>* pBoundaryCondition2 = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pBoundaryCondition2);
        
        // Linear solver
        SimpleNonlinearSolver solver;
        
        // Assembler
        SimpleNonlinearEllipticAssembler<1,1> assembler;
        
        // Set up solution guess for residuals
        int length=mesh.GetNumNodes();
                
        // Set up initial Guess
        Vec initialGuess;
        VecCreate(PETSC_COMM_WORLD, &initialGuess);
        VecSetSizes(initialGuess, PETSC_DECIDE,length);
        VecSetType(initialGuess, VECSEQ);
        for(int i=0; i<length ; i++)
        {
            VecSetValue(initialGuess, i, (0.5*((i/length)*(i/length)+1)), INSERT_VALUES);
        }
        VecAssemblyBegin(initialGuess);
        VecAssemblyEnd(initialGuess); 
        
        GaussianQuadratureRule<1> quadRule(2);
        LinearBasisFunction<1> basis_func;
        
        Vec answer;
        Vec residual;
        VecDuplicate(initialGuess,&residual);
        VecDuplicate(initialGuess,&answer);
        
//        Vec SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::AssembleSystem(
//                        ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh,
//                        AbstractNonlinearEllipticPde<SPACE_DIM> *pPde, 
//                        BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> *pBoundaryConditions,
//                        AbstractNonlinearSolver *pSolver,
//                        AbstractBasisFunction<SPACE_DIM> *pBasisFunction,
//                        GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule,
//                        Vec initialGuess,
//                        bool UseAnalyticalJacobian)

        try {
            answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, false);
        } catch (Exception e) {
            TS_TRACE(e.getMessage());
        }
        
        // Check result
        double *ans;
        int ierr = VecGetArray(answer, &ans);
        // Solution should be u = 0.5*x*(3-x)
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            double x = 0.0 + 0.01*i;
            double u = 0.5*(x*x+1);
            TS_ASSERT_DELTA(ans[i], u, 0.001);
        }
        VecRestoreArray(answer, &ans);
                
        // set up file output
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","CancerQ1");
        mpNewTestWriter->DefineFixedDimension("Space","dimensionless", mesh.GetNumElements()+1);
        int new_c_var_id = mpNewTestWriter->DefineVariable("Concentration","mM");
        int new_space_var_id = mpNewTestWriter->DefineVariable("Space","dimensionless");
        mpNewTestWriter->EndDefineMode();
        
        
        for (int i = 0; i < mesh.GetNumElements()+1; i++) 
        {
                 
            mpNewTestWriter->PutVariable(new_space_var_id, mesh.GetNodeAt(i)->GetPoint()[0], i);
            mpNewTestWriter->PutVariable(new_c_var_id, ans[i], i); 
        }
        mpNewTestWriter->Close();

    }
    
    void testQuestion2(void)
    {
        int FakeArgc=0;
        char *FakeArgv0="testrunner";
        char **FakeArgv=&FakeArgv0;
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
        
        // Create mesh from mesh reader
        
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        // Instantiate PDE object
        Practical1Question2Pde<1> pde;  
        
        // Boundary conditions
        // u'(0)=0, u(1)=1
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
        bcc.AddNeumannBoundaryCondition(*iter,pBoundaryCondition1);
        
        ConstBoundaryCondition<1>* pBoundaryCondition2 = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pBoundaryCondition2);
        
        // Linear solver
        SimpleNonlinearSolver solver;
        
        // Assembler
        SimpleNonlinearEllipticAssembler<1,1> assembler;
        
        // Set up solution guess for residuals
        int length=mesh.GetNumNodes();
                
        // Set up initial Guess
        Vec initialGuess;
        VecCreate(PETSC_COMM_WORLD, &initialGuess);
        VecSetSizes(initialGuess, PETSC_DECIDE,length);
        VecSetType(initialGuess, VECSEQ);
        for(int i=0; i<length ; i++)
        {
            VecSetValue(initialGuess, i, (0.5*((i/length)*(i/length)+1)), INSERT_VALUES);
        }
        VecAssemblyBegin(initialGuess);
        VecAssemblyEnd(initialGuess); 
        
        GaussianQuadratureRule<1> quadRule(2);
        LinearBasisFunction<1> basis_func;
        
        Vec answer;
        Vec residual;
        VecDuplicate(initialGuess,&residual);
        VecDuplicate(initialGuess,&answer);
        
        try {
            answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, &basis_func, &quadRule, initialGuess, true);
        } catch (Exception e) {
            TS_TRACE(e.getMessage());
        }
        
        // Check result
        double *ans;
        int ierr = VecGetArray(answer, &ans);
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(ans[i], 0.75, 0.250001);
            TS_ASSERT_DELTA(ans[i], 1.0, 0.250001);
        }
        VecRestoreArray(answer, &ans);
                
        // set up file output
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","CancerQ2a");
        mpNewTestWriter->DefineFixedDimension("Space","dimensionless", mesh.GetNumElements()+1);
        int new_c_var_id = mpNewTestWriter->DefineVariable("Concentration","mM");
        int new_space_var_id = mpNewTestWriter->DefineVariable("Space","dimensionless");
        mpNewTestWriter->EndDefineMode();
        
        
        for (int i = 0; i < mesh.GetNumElements()+1; i++) 
        {
                 
            mpNewTestWriter->PutVariable(new_space_var_id, mesh.GetNodeAt(i)->GetPoint()[0], i);
            mpNewTestWriter->PutVariable(new_c_var_id, ans[i], i); 
        }
        mpNewTestWriter->Close();

    }


    void testQuestions3to5(void)
    {
        int FakeArgc=0;
        char *FakeArgv0="testrunner";
        char **FakeArgv=&FakeArgv0;
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
        

                
         // Create mesh from mesh reader
        
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Solve for C
        
        // Instantiate PDE object
        Practical1Question2Pde<1> concPde;
        
        // Boundary conditions
        // u'(0)=0, u(1)=1
        BoundaryConditionsContainer<1,1> concBcc;
        ConstBoundaryCondition<1>* pConcBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator concIter = mesh.GetFirstBoundaryElement();
        concBcc.AddNeumannBoundaryCondition(*concIter,pConcBoundaryCondition1);
        
        ConstBoundaryCondition<1>* pConcBoundaryCondition2 = new ConstBoundaryCondition<1>(1.0);
        concBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pConcBoundaryCondition2);
        
        // Linear solver
        SimpleNonlinearSolver concSolver;
        
        // Assembler
        SimpleNonlinearEllipticAssembler<1,1> concAssembler;
        
        // Set up solution guess for residuals
        int length=mesh.GetNumNodes();
                
        // Set up initial Guess
        Vec initialGuess;
        VecCreate(PETSC_COMM_WORLD, &initialGuess);
        VecSetSizes(initialGuess, PETSC_DECIDE,length);
        VecSetType(initialGuess, VECSEQ);
        for(int i=0; i<length ; i++)
        {
            VecSetValue(initialGuess, i, (0.5*((i/length)*(i/length)+1)), INSERT_VALUES);
        }
        VecAssemblyBegin(initialGuess);
        VecAssemblyEnd(initialGuess); 
        
        GaussianQuadratureRule<1> quadRule(2);
        LinearBasisFunction<1> basis_func;
        
        Vec concAnswer;
        Vec concResidual;
        VecDuplicate(initialGuess,&concResidual);
        VecDuplicate(initialGuess,&concAnswer);
        
        try {
            concAnswer=concAssembler.AssembleSystem(&mesh, &concPde, &concBcc, &concSolver, &basis_func, &quadRule, initialGuess, false);
        } catch (Exception e) {
            TS_TRACE(e.getMessage());
        }
        
               
        // Check result
        double *concAns;
        int concIerr = VecGetArray(concAnswer, &concAns);
        // Solution should lie between 0.75 and 1.0
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(concAns[i], 0.75, 0.250001);
            TS_ASSERT_DELTA(concAns[i], 1.0, 0.250001);
        }
        VecRestoreArray(concAnswer, &concAns);
        
        // Calculate the value of x at which C = alpha
        double alpha = 0.8;
        int i = 0;
        double x_alpha = 0.0;
        while(concAns[i] < alpha && i<=mesh.GetNumElements()+1)
        {
            i++;
        }
        
        if (i > 0) 
        {
            x_alpha = mesh.GetNodeAt(0)->GetPoint()[0]  
                             + (mesh.GetNodeAt(100)->GetPoint()[0]-mesh.GetNodeAt(0)->GetPoint()[0])*(i-1)/( (double) length )
                             + (alpha - concAns[i-1])*(mesh.GetNodeAt(i)->GetPoint()[0]-mesh.GetNodeAt(i-1)->GetPoint()[0])/(concAns[i] - concAns[i-1]);
        }
   
        //std::cout << "x_alpha = " << x_alpha <<std::endl;
               
        
        // Solve for P_c
        
        // Instantiate PDE object
        Practical1Question3Pde<1> pressureCellPde;  
        
        pressureCellPde.mXalpha = x_alpha;
        
        // Boundary conditions
        // u'(0)=0, u(1)=1
        BoundaryConditionsContainer<1,1> pressureCellBcc;
        ConstBoundaryCondition<1>* pPressureCellBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator pressureCellIter = mesh.GetFirstBoundaryElement();
        pressureCellBcc.AddNeumannBoundaryCondition(*pressureCellIter,pPressureCellBoundaryCondition1);
        
        ConstBoundaryCondition<1>* pPressureCellBoundaryCondition2 = new ConstBoundaryCondition<1>(0.0);
        pressureCellBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pPressureCellBoundaryCondition2);
        
        // Linear solver
        SimpleLinearSolver pressureCellSolver;
        
        // Assembler
        SimpleLinearEllipticAssembler<1,1> pressureCellAssembler;
        
        Vec pressureCellAnswer = pressureCellAssembler.AssembleSystem(mesh, &pressureCellPde, pressureCellBcc, &pressureCellSolver);
        
        // Check result
        double *pressureCellAns;
        int pressureIerr = VecGetArray(pressureCellAnswer, &pressureCellAns);
        // Solution should lie between 1.2 and -1.0
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(pressureCellAns[i], -1.0, 2.20001);
            TS_ASSERT_DELTA(pressureCellAns[i], 1.2, 2.20001); 
        }
        VecRestoreArray(pressureCellAnswer, &pressureCellAns);
        
        
        // Solve for P_e
        
        double pressureExtAns[mesh.GetNumElements()+1];
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            pressureExtAns[i] = -pressureCellAns[i];
        }
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(pressureExtAns[i], -1.2, 2.20001);
            TS_ASSERT_DELTA(pressureExtAns[i], 1.0, 2.20001); 
        }
        
        // Solve for U_c
        
        double velocityCell[mesh.GetNumElements()+1];
        
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            if (i == 0)
            {
                // Use forward difference
                velocityCell[i] = - (pressureCellAns[i+1] - pressureCellAns[i])*mesh.GetNumElements();
            }
            else if(i == mesh.GetNumElements()+1)
            {
                // Use backward difference
                velocityCell[i] = - (pressureCellAns[i] - pressureCellAns[i-1])*mesh.GetNumElements();
            }
            else
            {
                // Use central difference
                velocityCell[i] = - (pressureCellAns[i+1] - pressureCellAns[i-1])*mesh.GetNumElements()/2.0;
            }
            
        }
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(pressureCellAns[i], -0.4, 0.70001);
            TS_ASSERT_DELTA(pressureCellAns[i], 0.3, 0.70001); 
        }
              
        // Solve for U_e
        double phi0 = 0.5;
        double velocityExt[mesh.GetNumElements()+1];
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            velocityExt[i] = (-phi0/(1-phi0))*velocityCell[i];
        }
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(pressureCellAns[i], -0.3, 0.70001);
            TS_ASSERT_DELTA(pressureCellAns[i], 0.4, 0.70001); 
        }
                
        BoundaryOdeSystem* pBoundarySystem = new BoundaryOdeSystem();
        EulerIvpOdeSolver* myEulerSolver = new EulerIvpOdeSolver;
        OdeSolution solutions;
        
        std::vector<double> boundaryInit(1);

        boundaryInit[0] = 1.0;
         
        pBoundarySystem->mVelocityCellAtBoundary = velocityCell[mesh.GetNumElements()+1];
        
        solutions = myEulerSolver->Solve(pBoundarySystem, 0.0, 0.01, 0.01, boundaryInit);  
        
        
        // set up file output
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","CancerQ5");
        mpNewTestWriter->DefineFixedDimension("Space","dimensionless", mesh.GetNumElements()+1);
        int new_c_var_id = mpNewTestWriter->DefineVariable("Concentration","mM");
        int new_space_var_id = mpNewTestWriter->DefineVariable("Space","dimensionless");
        int new_x_var_id = mpNewTestWriter->DefineVariable("X","dimensionless");
        int new_p_var_id = mpNewTestWriter->DefineVariable("Pressure","Pa");
        mpNewTestWriter->EndDefineMode();
        
        
        for (int i = 0; i < mesh.GetNumElements()+1; i++) 
        {
                 
            mpNewTestWriter->PutVariable(new_space_var_id, mesh.GetNodeAt(i)->GetPoint()[0], i);
            mpNewTestWriter->PutVariable(new_p_var_id, pressureCellAns[i], i);
            mpNewTestWriter->PutVariable(new_x_var_id, mesh.GetNodeAt(100)->GetPoint()[0], i);
            mpNewTestWriter->PutVariable(new_c_var_id, concAns[i], i); 
        }
        mpNewTestWriter->Close();
   
   
        
    }

};

#endif //_TESTPRACTICALONE_HPP_
