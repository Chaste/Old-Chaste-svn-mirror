#ifndef _TESTPRACTICALONE_HPP_
#define _TESTPRACTICALONE_HPP_

#include "SimpleLinearSolver.hpp"
#include "petscmat.h"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "SimpleNonlinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.hpp"
#include <vector>
#include <iostream>
#include "Node.hpp" 
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ColumnDataWriter.hpp"

#include "Practical1Question1Pde.hpp"
#include "Practical1Question1PdeNonlinear.hpp"
#include "ConcentrationPde.hpp"
#include "CellPressurePde.hpp"
#include "Practical1Question2FudgePde.hpp"
#include "Practical1Question3FudgePde.hpp"
#include "BoundaryOdeSystem.hpp"

#include "SimpleLinearEllipticAssembler.hpp"
#include "SimpleNonlinearEllipticAssembler.hpp"

#include "EulerIvpOdeSolver.hpp"

  
class TestPracticalOne : public CxxTest::TestSuite 
{
public:
	
	void testPrac1Question1(void)
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
    
    void testPrac1Question1Nonlinear(void)
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
                
        Vec answer;
        Vec residual;
        VecDuplicate(initialGuess,&residual);
        VecDuplicate(initialGuess,&answer);

        try {
            answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initialGuess, false);
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
                
//        // set up file output
//        ColumnDataWriter *mpNewTestWriter;
//        mpNewTestWriter = new ColumnDataWriter("data","CancerQ1");
//        mpNewTestWriter->DefineFixedDimension("Space","dimensionless", mesh.GetNumElements()+1);
//        int new_c_var_id = mpNewTestWriter->DefineVariable("Concentration","mM");
//        int new_space_var_id = mpNewTestWriter->DefineVariable("Space","dimensionless");
//        mpNewTestWriter->EndDefineMode();
//        
//        
//        for (int i = 0; i < mesh.GetNumElements()+1; i++) 
//        {
//                 
//            mpNewTestWriter->PutVariable(new_space_var_id, mesh.GetNodeAt(i)->GetPoint()[0], i);
//            mpNewTestWriter->PutVariable(new_c_var_id, ans[i], i); 
//        }
//        mpNewTestWriter->Close();

    }
    
    void testPrac1Question2(void)
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
        ConcentrationPde<1> pde;  
        
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
        
        Vec answer;
        Vec residual;
        VecDuplicate(initialGuess,&residual);
        VecDuplicate(initialGuess,&answer);
        
        try {
            answer=assembler.AssembleSystem(&mesh, &pde, &bcc, &solver, initialGuess, true);
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
                
//        // set up file output
//        ColumnDataWriter *mpNewTestWriter;
//        mpNewTestWriter = new ColumnDataWriter("data","CancerQ2a");
//        mpNewTestWriter->DefineFixedDimension("Space","dimensionless", mesh.GetNumElements()+1);
//        int new_c_var_id = mpNewTestWriter->DefineVariable("Concentration","mM");
//        int new_space_var_id = mpNewTestWriter->DefineVariable("Space","dimensionless");
//        mpNewTestWriter->EndDefineMode();
//        
//        
//        for (int i = 0; i < mesh.GetNumElements()+1; i++) 
//        {
//                 
//            mpNewTestWriter->PutVariable(new_space_var_id, mesh.GetNodeAt(i)->GetPoint()[0], i);
//            mpNewTestWriter->PutVariable(new_c_var_id, ans[i], i); 
//        }
//        mpNewTestWriter->Close();

    }


    void testPrac1Questions3to5(void)
    {
        int FakeArgc=0;
        char *FakeArgv0="testrunner";
        char **FakeArgv=&FakeArgv0;
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
        
        double alpha = 0.8;
                
         // Create mesh from mesh reader
        
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        
        // Instantiate PDE object for oxygen concentration
        ConcentrationPde<1> concentrationPde;
        // Instantiate PDE object for cell Pressure
        CellPressurePde<1> cellPressurePde;  
        
        // Boundary conditions for oxygen concentration
        // C'(0)=0, C(1)=1
        BoundaryConditionsContainer<1,1> concentrationBcc;
        ConstBoundaryCondition<1>* pConcentrationBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator concentrationIter = mesh.GetFirstBoundaryElement();
        concentrationBcc.AddNeumannBoundaryCondition(*concentrationIter,pConcentrationBoundaryCondition1);
        
        ConstBoundaryCondition<1>* pConcentrationBoundaryCondition2 = new ConstBoundaryCondition<1>(1.0);
        concentrationBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pConcentrationBoundaryCondition2);
        
        // Boundary conditions for cell pressure
        // Pc'(0)=0, Pc(1)=0
        BoundaryConditionsContainer<1,1> cellPressureBcc;
        ConstBoundaryCondition<1>* pCellPressureBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator cellPressureIter = mesh.GetFirstBoundaryElement();
        cellPressureBcc.AddNeumannBoundaryCondition(*cellPressureIter,pCellPressureBoundaryCondition1);
        
        ConstBoundaryCondition<1>* pCellPressureBoundaryCondition2 = new ConstBoundaryCondition<1>(0.0);
        cellPressureBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pCellPressureBoundaryCondition2);
        
        //Instantiate ODE system for external tumour boundary
        BoundaryOdeSystem* pBoundarySystem = new BoundaryOdeSystem();
        EulerIvpOdeSolver* myEulerSolver = new EulerIvpOdeSolver;
        OdeSolution solutions;
        // Boundary Conditions for ODE on boundary
        std::vector<double> boundaryInit(1);

        boundaryInit[0] = 1.0;
        
        
        // Nonlinear solver for oxygen concentration
        SimpleNonlinearSolver concentrationSolver;
        
        //Linear Solver for Cell Pressure
        SimpleLinearSolver cellPressureSolver;
        
        
        // Assembler for oxygen concentration
        SimpleNonlinearEllipticAssembler<1,1> concentrationAssembler;
        
        // Assembler for cell pressure
        SimpleLinearEllipticAssembler<1,1> cellPressureAssembler;
        
        
        
        
        
        
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
        
        Vec concentrationAnswer;
        Vec concentrationResidual;
        VecDuplicate(initialGuess,&concentrationResidual);
        VecDuplicate(initialGuess,&concentrationAnswer);
        
        try {
            concentrationAnswer=concentrationAssembler.AssembleSystem(&mesh, &concentrationPde, &concentrationBcc, &concentrationSolver, initialGuess, false);
        } catch (Exception e) {
            TS_TRACE(e.getMessage());
        }
        
               
        // Check result
        double *concentrationAns;
        int concentrationIerr = VecGetArray(concentrationAnswer, &concentrationAns);
        // Solution should lie between 0.75 and 1.0
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(concentrationAns[i], 0.75, 0.250001);
            TS_ASSERT_DELTA(concentrationAns[i], 1.0, 0.250001);
        }
        VecRestoreArray(concentrationAnswer, &concentrationAns);
        
        // Calculate the value of x at which C = alpha , this boundary determines where the cells begin to die of oxygen deprivation
        
        int i = 0;
        double x_alpha = 0.0;
        while(concentrationAns[i] < alpha && i<=mesh.GetNumElements()+1)
        {
            i++;
        }
        
        if (i > 0) 
        {
            x_alpha = mesh.GetNodeAt(0)->GetPoint()[0]  
                             + (mesh.GetNodeAt(100)->GetPoint()[0]-mesh.GetNodeAt(0)->GetPoint()[0])*(i-1)/( (double) length )
                             + (alpha - concentrationAns[i-1])*(mesh.GetNodeAt(i)->GetPoint()[0]-mesh.GetNodeAt(i-1)->GetPoint()[0])/(concentrationAns[i] - concentrationAns[i-1]);
        }
         cellPressurePde.mXalpha = x_alpha;
        
        // Solve the linear pde for cell Pressure                
        Vec cellPressureAnswer = cellPressureAssembler.AssembleSystem(mesh, &cellPressurePde, cellPressureBcc, &cellPressureSolver);
        
        // Test to see if result lies within bounds calculated from Matlab
        double *cellPressureAns;
        int cellPressureIerr = VecGetArray(cellPressureAnswer, &cellPressureAns);
        // Solution should lie between 0.06 and -0.07
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(cellPressureAns[i], 0.07, 0.14);
            TS_ASSERT_DELTA(cellPressureAns[i], -0.07, 0.14); 
        }
        VecRestoreArray(cellPressureAnswer, &cellPressureAns);
        
        
        // Solve for extracellular fluid pressure
        
        double fluidPressureAns[mesh.GetNumElements()+1];
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            fluidPressureAns[i] = -cellPressureAns[i];
        }
        
        // Solution should lie between -0.06 and 0.07 ,bounds obtained from matlab
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(fluidPressureAns[i], -0.07, 0.14);
            TS_ASSERT_DELTA(fluidPressureAns[i], 0.07, 0.14); 
        }
        
        // Solve for cell velocity , U_c partial derivative of cell pressure w.r.t x
        
        double cellVelocity[mesh.GetNumElements()+1];
        
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            if (i == 0)
            {
                // Use forward difference
                cellVelocity[i] = - (cellPressureAns[i+1] - cellPressureAns[i])*mesh.GetNumElements();
            }
            else if(i == mesh.GetNumElements())
            {
                // Use backward difference
                cellVelocity[i] = - (cellPressureAns[i] - cellPressureAns[i-1])*mesh.GetNumElements();
            }
            else
            {
                // Use central difference
                cellVelocity[i] = - (cellPressureAns[i+1] - cellPressureAns[i-1])*mesh.GetNumElements()/2.0;
            }
            
        }
        
        //Test cell Pressure lies within bounds from matlab
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(cellPressureAns[i], -0.4, 0.70001);
            TS_ASSERT_DELTA(cellPressureAns[i], 0.3, 0.70001); 
        }
              
        // Solve for U_e
        double phi0 = 0.5;
        double fluidVelocity[mesh.GetNumElements()+1];
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            fluidVelocity[i] = (-phi0/(1-phi0))*cellVelocity[i];
        }
        for (int i=0; i < mesh.GetNumElements()+1; i++)
        {
            TS_ASSERT_DELTA(cellPressureAns[i], -0.3, 0.70001);
            TS_ASSERT_DELTA(cellPressureAns[i], 0.4, 0.70001); 
        }
                
        
         
        pBoundarySystem->mCellVelocityAtBoundary = cellVelocity[mesh.GetNumElements()+1];
        
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
            mpNewTestWriter->PutVariable(new_p_var_id, cellPressureAns[i], i);
            mpNewTestWriter->PutVariable(new_x_var_id, mesh.GetNodeAt(100)->GetPoint()[0], i);
            mpNewTestWriter->PutVariable(new_c_var_id, concentrationAns[i], i); 
        }
        mpNewTestWriter->Close();
   
        //read in good data file and compare line by line
        std::ifstream testfile("data/CancerQ5.dat",std::ios::in);
        std::ifstream goodfile("data/CancerQ5Good.dat",std::ios::in);
        std::string teststring;
        std::string goodstring;
        while(getline(testfile, teststring))
        {
              getline(goodfile,goodstring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();
        
    }


    void testPrac1Question6(void)
    {
        
        ///////////////////////////////////////////////
        // Set up 
        ///////////////////////////////////////////////
        
        int FakeArgc=0;
        char *FakeArgv0="testrunner";
        char **FakeArgv=&FakeArgv0;
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);

        // We consider only the 1-D case
        const int space_dim = 1;
        const int element_dim = 1;
        
        double time = 0.0;
        double timestep = 0.01;
        int endTimeIndex = 100;         // Number of timesteps
        
        // Initialise tumour boundary
        double X_new = 1.0;
        
        // Define alpha and initialise x_alpha 
        const double alpha = 0.5;
        double x_alpha = 0.0;
        
        double temp;
                
        // Initialise the scale factor        
        double scaleFactor = 0.0;
        
                
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
        ConformingTetrahedralMesh<element_dim,space_dim> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object for oxygen concentration
        ConcentrationPde<space_dim> concentrationPde;
        // Instantiate PDE object for cell Pressure
        CellPressurePde<space_dim> cellPressurePde;  
        
        //Instantiate ODE system for external tumour boundary
        BoundaryOdeSystem* pBoundarySystem = new BoundaryOdeSystem();
        EulerIvpOdeSolver* myEulerSolver = new EulerIvpOdeSolver;
        OdeSolution solutions;
        
        // Boundary Conditions for ODE for tumour boundary
        std::vector<double> boundaryInit(1);
        boundaryInit[0] = X_new;
        
        // Nonlinear solver for oxygen concentration
        SimpleNonlinearSolver concentrationSolver;
        
        //Linear Solver for Cell Pressure
        SimpleLinearSolver cellPressureSolver;
        
        // Assembler for oxygen concentration
        SimpleNonlinearEllipticAssembler<element_dim,space_dim> concentrationAssembler;
        
        // Assembler for cell pressure
        SimpleLinearEllipticAssembler<element_dim,space_dim> cellPressureAssembler;
        
        // Set up solution guess for residuals
        int length=mesh.GetNumNodes();
                
        // Set up initial guess for concentration
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
        
        Vec concentrationAnswer;
        Vec concentrationResidual;
        VecDuplicate(initialGuess,&concentrationResidual);
        VecDuplicate(initialGuess,&concentrationAnswer);
        
        // Array to hold cell pressure solution
        double *cellPressureAns;
        
        // Set up data writer
        ColumnDataWriter *mpTimeWriter;
        mpTimeWriter = new ColumnDataWriter("data","CancerQ6Time");
        mpTimeWriter->DefineUnlimitedDimension("Time","s");
        int time_id = mpTimeWriter->DefineVariable("Time","dim.");
        int xt_var_id = mpTimeWriter->DefineVariable("X(t)","dim.");
        mpTimeWriter->EndDefineMode();
        
        
        ///////////////////////////////////////////////               
        // Start of grand time loop
        ///////////////////////////////////////////////
        
        for(int timeIndex = 0 ; timeIndex < endTimeIndex ; timeIndex++)
        {
            
        // Boundary conditions for oxygen concentration
        // C'(0)=0, C(1)=1
        BoundaryConditionsContainer<element_dim,space_dim> concentrationBcc;
        ConstBoundaryCondition<space_dim>* pConcentrationBoundaryCondition1 = new ConstBoundaryCondition<space_dim>(0.0);
        
        ConformingTetrahedralMesh<element_dim,space_dim>::BoundaryElementIterator concentrationIter = mesh.GetFirstBoundaryElement();
        concentrationBcc.AddNeumannBoundaryCondition(*concentrationIter,pConcentrationBoundaryCondition1);
        
        ConstBoundaryCondition<space_dim>* pConcentrationBoundaryCondition2 = new ConstBoundaryCondition<space_dim>(1.0);
        concentrationBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pConcentrationBoundaryCondition2);
        
        // Boundary conditions for cell pressure
        // Pc'(0)=0, Pc(1)=0
        BoundaryConditionsContainer<element_dim,space_dim> cellPressureBcc;
        ConstBoundaryCondition<space_dim>* pCellPressureBoundaryCondition1 = new ConstBoundaryCondition<space_dim>(0.0);
        
        ConformingTetrahedralMesh<element_dim,space_dim>::BoundaryElementIterator cellPressureIter = mesh.GetFirstBoundaryElement();
        cellPressureBcc.AddNeumannBoundaryCondition(*cellPressureIter,pCellPressureBoundaryCondition1);
        
        ConstBoundaryCondition<space_dim>* pCellPressureBoundaryCondition2 = new ConstBoundaryCondition<space_dim>(0.0);
        cellPressureBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pCellPressureBoundaryCondition2);
        
            
            try {
                concentrationAnswer=concentrationAssembler.AssembleSystem(&mesh, &concentrationPde, &concentrationBcc, &concentrationSolver, initialGuess, true);
            } catch (Exception e) {
                TS_TRACE(e.getMessage());
            }
                               
            // Check result
            double *concentrationAns;
            int concentrationIerr = VecGetArray(concentrationAnswer, &concentrationAns);
            
            // The Solver fails if the initial guess is updated!!
            
            // Calculate the value of x at which C = alpha , this boundary determines where the cells begin to die of oxygen deprivation
            
            int i = 0;
            
            while(concentrationAns[i] < alpha && i<=mesh.GetNumElements()+1)
            {
                i++;
            }
            
            if (i > 0) 
            {
                x_alpha = mesh.GetNodeAt(0)->GetPoint()[0]  
                                 + (mesh.GetNodeAt(100)->GetPoint()[0]-mesh.GetNodeAt(0)->GetPoint()[0])*(i-1)/( (double) length )
                                 + (alpha - concentrationAns[i-1])*(mesh.GetNodeAt(i)->GetPoint()[0]-mesh.GetNodeAt(i-1)->GetPoint()[0])/(concentrationAns[i] - concentrationAns[i-1]);
            }
             cellPressurePde.mXalpha = x_alpha;
            //std::cout<<"x_alpha = " << x_alpha << std::endl;
            // Solve the linear pde for cell Pressure                
            Vec cellPressureAnswer = cellPressureAssembler.AssembleSystem(mesh, &cellPressurePde, cellPressureBcc, &cellPressureSolver);
            
//            // Test to see if result lies within bounds calculated from Matlab
            
            int cellPressureIerr = VecGetArray(cellPressureAnswer, &cellPressureAns);
//            // Solution should lie between 0.06 and -0.07
//            for (int i=0; i < mesh.GetNumElements()+1; i++)
//            {
//                TS_ASSERT_DELTA(cellPressureAns[i], 0.07, 0.14);
//                TS_ASSERT_DELTA(cellPressureAns[i], -0.07, 0.14); 
//            }
            VecRestoreArray(cellPressureAnswer, &cellPressureAns);
            
            //std::cout<< " cell Pressure " << cellPressureAns[99] << "  \n"  ;
            // Solve for extracellular fluid pressure
            
            double fluidPressureAns[mesh.GetNumElements()+1];
            for (int i=0; i < mesh.GetNumElements()+1; i++)
            {
                fluidPressureAns[i] = -cellPressureAns[i];
            }
            
            
//            // Solution should lie between -0.06 and 0.07 ,bounds obtained from matlab
//            for (int i=0; i < mesh.GetNumElements()+1; i++)
//            {
//                TS_ASSERT_DELTA(fluidPressureAns[i], -0.07, 0.14);
//                TS_ASSERT_DELTA(fluidPressureAns[i], 0.07, 0.14); 
//            }
            
            // Solve for cell velocity , U_c partial derivative of cell pressure w.r.t x
            
            double cellVelocity[mesh.GetNumElements()+1];
            
                       
            for (int i=0; i < mesh.GetNumElements()+1; i++)
            {
                if (i == 0)
                {
                    // Use forward difference
                    cellVelocity[i] = - (cellPressureAns[i+1] - cellPressureAns[i])*mesh.GetNumElements()/X_new;
                }
                else if(i == mesh.GetNumElements())
                {
                    // Use backward difference
                    cellVelocity[i] = - (cellPressureAns[i] - cellPressureAns[i-1])*mesh.GetNumElements()/X_new;
                }
                else
                {
                    // Use central difference
                    cellVelocity[i] = - (cellPressureAns[i+1] - cellPressureAns[i-1])*mesh.GetNumElements()/(2.0*X_new);
                }
                
            }
            
//            //Test cell pressure lies within bounds from matlab
//            for (int i=0; i < mesh.GetNumElements()+1; i++)
//            {
//                TS_ASSERT_DELTA(cellPressureAns[i], -0.4, 0.70001);
//                TS_ASSERT_DELTA(cellPressureAns[i], 0.3, 0.70001); 
//            }
                  
            // Solve for U_e
            double phi0 = 0.5;
            double fluidVelocity[mesh.GetNumElements()+1];
            for (int i=0; i < mesh.GetNumElements()+1; i++)
            {
                fluidVelocity[i] = (-phi0/(1-phi0))*cellVelocity[i];
            }
//            for (int i=0; i < mesh.GetNumElements()+1; i++)
//            {
//                TS_ASSERT_DELTA(cellPressureAns[i], -0.3, 0.70001);
//                TS_ASSERT_DELTA(cellPressureAns[i], 0.4, 0.70001); 
//            }

            //std::cout << "Cell Velocity at boundary = " <<cellVelocity[mesh.GetNumElements()]<< "\n";                         
            pBoundarySystem->mCellVelocityAtBoundary = cellVelocity[mesh.GetNumElements()];
            boundaryInit[0] = X_new;
            solutions = myEulerSolver->Solve(pBoundarySystem, time, time+timestep, timestep, boundaryInit);  
            
            // Rescale Mesh
            // Get X_new, rescale with X_old
            X_new = solutions.mSolutions[1][0];
            //std::cout << "X_new = "<< X_new << std::endl;
            scaleFactor = solutions.mSolutions[1][0] / solutions.mSolutions[0][0];
            //std::cout << "Scale Factor = " << scaleFactor << std::endl;
            Point<1> newBoundaryLocation(X_new);
            
            mesh.RescaleMeshFromBoundaryNode(newBoundaryLocation, mesh.GetNumElements());
            
            //std::cout << "New RHS node location = " << mesh.GetNodeAt(mesh.GetNumNodes()-1)->GetPoint()[0] << std::endl;           
            time = time + timestep;
        // End of grand time loop
            if(timeIndex==23)
            {
                    // set up file output
                ColumnDataWriter *mpNewTestWriter;
                mpNewTestWriter = new ColumnDataWriter("data","CancerQ6t23");
                mpNewTestWriter->DefineFixedDimension("Space","dimensionless", mesh.GetNumElements()+1);
                int new_c_var_id = mpNewTestWriter->DefineVariable("Concentration","mM");
                int new_space_var_id = mpNewTestWriter->DefineVariable("Space","dimensionless");
                int new_x_var_id = mpNewTestWriter->DefineVariable("X","dimensionless");
                int new_p_var_id = mpNewTestWriter->DefineVariable("Pressure","Pa");
                int uc_var_id = mpNewTestWriter->DefineVariable("Cell Velocity"," ");
                mpNewTestWriter->EndDefineMode();
                
                for (int i = 0; i < mesh.GetNumElements()+1; i++) 
                {
                         
                    mpNewTestWriter->PutVariable(new_space_var_id, mesh.GetNodeAt(i)->GetPoint()[0], i);
                    mpNewTestWriter->PutVariable(new_p_var_id, cellPressureAns[i], i);
                    mpNewTestWriter->PutVariable(new_x_var_id, mesh.GetNodeAt(100)->GetPoint()[0], i);
                    mpNewTestWriter->PutVariable(new_c_var_id, concentrationAns[i], i); 
                    mpNewTestWriter->PutVariable(uc_var_id, cellVelocity[i], i); 
                }
                mpNewTestWriter->Close();
           
            }
            mpTimeWriter->PutVariable(time_id, time);
            mpTimeWriter->PutVariable(xt_var_id, X_new);
            mpTimeWriter->AdvanceAlongUnlimitedDimension();     
                
        }

        mpTimeWriter->Close();
    }
//
//
//    void donttestQuestion6MeshRescaleFudge(void)
//    {
//        int FakeArgc=0;
//        char *FakeArgv0="testrunner";
//        char **FakeArgv=&FakeArgv0;
//        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
//        
//        const int space_dim = 1;
//        const int element_dim = 1;
//        double alpha = 0.8;
//        int endTimeIndex = 20;
//        double time = 0.0;
//        double timestep = 0.01;
//        double temp;
//        double X_new = 1.0;
//                
//         // Create mesh from mesh reader
//        
//        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
//        ConformingTetrahedralMesh<element_dim,space_dim> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//        
//           
//        // Instantiate PDE object for oxygen concentration
//        Practical1Question2FudgePde<space_dim> concentrationPde;
//        // Instantiate PDE object for cell Pressure
//        Practical1Question3FudgePde<space_dim> cellPressurePde;
//        
//        concentrationPde.mXnew = X_new;
//        cellPressurePde.mXnew = X_new;
//        
//        // Boundary conditions for oxygen concentration
//        // C'(0)=0, C(1)=1
//        BoundaryConditionsContainer<element_dim,space_dim> concentrationBcc;
//        ConstBoundaryCondition<space_dim>* pConcentrationBoundaryCondition1 = new ConstBoundaryCondition<space_dim>(0.0);
//        
//        ConformingTetrahedralMesh<element_dim,space_dim>::BoundaryElementIterator concentrationIter = mesh.GetFirstBoundaryElement();
//        concentrationBcc.AddNeumannBoundaryCondition(*concentrationIter,pConcentrationBoundaryCondition1);
//        
//        ConstBoundaryCondition<space_dim>* pConcentrationBoundaryCondition2 = new ConstBoundaryCondition<space_dim>(1.0);
//        concentrationBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pConcentrationBoundaryCondition2);
//        
//        // Boundary conditions for cell pressure
//        // Pc'(0)=0, Pc(1)=0
//        BoundaryConditionsContainer<element_dim,space_dim> cellPressureBcc;
//        ConstBoundaryCondition<space_dim>* pCellPressureBoundaryCondition1 = new ConstBoundaryCondition<space_dim>(0.0);
//        
//        ConformingTetrahedralMesh<element_dim,space_dim>::BoundaryElementIterator cellPressureIter = mesh.GetFirstBoundaryElement();
//        cellPressureBcc.AddNeumannBoundaryCondition(*cellPressureIter,pCellPressureBoundaryCondition1);
//        
//        ConstBoundaryCondition<space_dim>* pCellPressureBoundaryCondition2 = new ConstBoundaryCondition<space_dim>(0.0);
//        cellPressureBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pCellPressureBoundaryCondition2);
//        
//        //Instantiate ODE system for external tumour boundary
//        BoundaryOdeSystem* pBoundarySystem = new BoundaryOdeSystem();
//        EulerIvpOdeSolver* myEulerSolver = new EulerIvpOdeSolver;
//        OdeSolution solutions;
//        // Boundary Conditions for ODE on boundary
//        std::vector<double> boundaryInit(1);
//
//        boundaryInit[0] = 1.0;
//        
//        
//        // Nonlinear solver for oxygen concentration
//        SimpleNonlinearSolver concentrationSolver;
//        
//        //Linear Solver for Cell Pressure
//        SimpleLinearSolver cellPressureSolver;
//        
//        
//        // Assembler for oxygen concentration
//        SimpleNonlinearEllipticAssembler<element_dim,space_dim> concentrationAssembler;
//        
//        // Assembler for cell pressure
//        SimpleLinearEllipticAssembler<element_dim,space_dim> cellPressureAssembler;
//        
//                
//        // Set up solution guess for residuals
//        int length=mesh.GetNumNodes();
//                
//        // Set up initial Guess
//        Vec initialGuess;
//        VecCreate(PETSC_COMM_WORLD, &initialGuess);
//        VecSetSizes(initialGuess, PETSC_DECIDE,length);
//        VecSetType(initialGuess, VECSEQ);
//        for(int i=0; i<length ; i++)
//        {
//            VecSetValue(initialGuess, i, (0.5*((i/length)*(i/length)+1)), INSERT_VALUES);
//        }
//        VecAssemblyBegin(initialGuess);
//        VecAssemblyEnd(initialGuess); 
//        
//        GaussianQuadratureRule<element_dim> quadRule(2);
//        LinearBasisFunction<element_dim> basis_func;
//        
//        Vec concentrationAnswer;
//        Vec concentrationResidual;
//        VecDuplicate(initialGuess,&concentrationResidual);
//        VecDuplicate(initialGuess,&concentrationAnswer);
//        
//        double scaleFactor = 0.0;
//        // Start of grand time loop
//        for(int timeIndex = 0 ; timeIndex < endTimeIndex ; timeIndex++)
//        {
//            
//            try {
//                concentrationAnswer=concentrationAssembler.AssembleSystem(&mesh, &concentrationPde, &concentrationBcc, &concentrationSolver, &basis_func, &quadRule, initialGuess, true);
//            } catch (Exception e) {
//                TS_TRACE(e.getMessage());
//            }
//                               
//            // Check result
//            double *concentrationAns;
//            int concentrationIerr = VecGetArray(concentrationAnswer, &concentrationAns);
//            
//            // The Solver fails if the initial guess is updated!!
////            // Update initial guess for next time
////            for(int i=0; i<length ; i++)
////            {
////                VecSetValue(initialGuess, i, concentrationAns[i], INSERT_VALUES);
////            }
////            VecAssemblyBegin(initialGuess);
////            VecAssemblyEnd(initialGuess); 
//            
//            // Calculate the value of x at which C = alpha , this boundary determines where the cells begin to die of oxygen deprivation
//            
//            int i = 0;
//            double x_alpha = 0.0;
//            while(concentrationAns[i] < alpha && i<=mesh.GetNumElements()+1)
//            {
//                i++;
//            }
//            
//            if (i > 0) 
//            {
//                x_alpha = mesh.GetNodeAt(0)->GetPoint()[0]  
//                                 + (mesh.GetNodeAt(100)->GetPoint()[0]-mesh.GetNodeAt(0)->GetPoint()[0])*(i-1)/( (double) length )
//                                 + (alpha - concentrationAns[i-1])*(mesh.GetNodeAt(i)->GetPoint()[0]-mesh.GetNodeAt(i-1)->GetPoint()[0])/(concentrationAns[i] - concentrationAns[i-1]);
//            }
//            cellPressurePde.mXalpha = x_alpha/X_new;
//            
//            // Solve the linear pde for cell Pressure                
//            Vec cellPressureAnswer = cellPressureAssembler.AssembleSystem(mesh, &cellPressurePde, cellPressureBcc, &cellPressureSolver);
//            
////            // Test to see if result lies within bounds calculated from Matlab
//            double *cellPressureAns;
//            int cellPressureIerr = VecGetArray(cellPressureAnswer, &cellPressureAns);
////            // Solution should lie between 0.06 and -0.07
////            for (int i=0; i < mesh.GetNumElements()+1; i++)
////            {
////                TS_ASSERT_DELTA(cellPressureAns[i], 0.07, 0.14);
////                TS_ASSERT_DELTA(cellPressureAns[i], -0.07, 0.14); 
////            }
//            VecRestoreArray(cellPressureAnswer, &cellPressureAns);
//            
//            
//            // Solve for extracellular fluid pressure
//            
//            double fluidPressureAns[mesh.GetNumElements()+1];
//            for (int i=0; i < mesh.GetNumElements()+1; i++)
//            {
//                fluidPressureAns[i] = -cellPressureAns[i];
//            }
//            
//            
////            // Solution should lie between -0.06 and 0.07 ,bounds obtained from matlab
////            for (int i=0; i < mesh.GetNumElements()+1; i++)
////            {
////                TS_ASSERT_DELTA(fluidPressureAns[i], -0.07, 0.14);
////                TS_ASSERT_DELTA(fluidPressureAns[i], 0.07, 0.14); 
////            }
//            
//            // Solve for cell velocity , U_c partial derivative of cell pressure w.r.t x
//            
//            double cellVelocity[mesh.GetNumElements()+1];
//            
//                       
//            for (int i=0; i < mesh.GetNumElements()+1; i++)
//            {
//                if (i == 0)
//                {
//                    // Use forward difference
//                    cellVelocity[i] = - X_new*(cellPressureAns[i+1] - cellPressureAns[i])*mesh.GetNumElements();
//                }
//                else if(i == mesh.GetNumElements())
//                {
//                    // Use backward difference
//                    cellVelocity[i] = - X_new*(cellPressureAns[i] - cellPressureAns[i-1])*mesh.GetNumElements();
//                }
//                else
//                {
//                    // Use central difference
//                    cellVelocity[i] = - X_new*(cellPressureAns[i+1] - cellPressureAns[i-1])*mesh.GetNumElements()/2.0;
//                }
//                
//            }
//            
////            //Test cell pressure lies within bounds from matlab
////            for (int i=0; i < mesh.GetNumElements()+1; i++)
////            {
////                TS_ASSERT_DELTA(cellPressureAns[i], -0.4, 0.70001);
////                TS_ASSERT_DELTA(cellPressureAns[i], 0.3, 0.70001); 
////            }
//                  
//            // Solve for U_e
//            double phi0 = 0.5;
//            double fluidVelocity[mesh.GetNumElements()+1];
//            for (int i=0; i < mesh.GetNumElements()+1; i++)
//            {
//                fluidVelocity[i] = (-phi0/(1-phi0))*cellVelocity[i];
//            }
////            for (int i=0; i < mesh.GetNumElements()+1; i++)
////            {
////                TS_ASSERT_DELTA(cellPressureAns[i], -0.3, 0.70001);
////                TS_ASSERT_DELTA(cellPressureAns[i], 0.4, 0.70001); 
////            }
//
//            //std::cout << "Cell Velocity = " <<cellVelocity[mesh.GetNumElements()]<< "\n";                         
//            pBoundarySystem->mCellVelocityAtBoundary = cellVelocity[mesh.GetNumElements()];
//            
//            solutions = myEulerSolver->Solve(pBoundarySystem, time, time+timestep, timestep, boundaryInit);  
//            
//            // Rescale Mesh
//            // Get X_new, rescale with X_old
//            
//            
//            X_new = solutions.mSolutions[1][0];
//            concentrationPde.mXnew = X_new;
//            cellPressurePde.mXnew = X_new;
//            std::cout << "X_new = "<< X_new << std::endl;
//            
//            //std::cout << "New RHS node location = " << mesh.GetNodeAt(mesh.GetNumNodes()-1)->GetPoint()[0] << std::endl;           
//            time = time + timestep;
//        // End of grand time loop
//        }    
//        
////        // set up file output
////        ColumnDataWriter *mpNewTestWriter;
////        mpNewTestWriter = new ColumnDataWriter("data","CancerQ5");
////        mpNewTestWriter->DefineFixedDimension("Space","dimensionless", mesh.GetNumElements()+1);
////        int new_c_var_id = mpNewTestWriter->DefineVariable("Concentration","mM");
////        int new_space_var_id = mpNewTestWriter->DefineVariable("Space","dimensionless");
////        int new_x_var_id = mpNewTestWriter->DefineVariable("X","dimensionless");
////        int new_p_var_id = mpNewTestWriter->DefineVariable("Pressure","Pa");
////        mpNewTestWriter->EndDefineMode();
////        
////        
////        for (int i = 0; i < mesh.GetNumElements()+1; i++) 
////        {
////                 
////            mpNewTestWriter->PutVariable(new_space_var_id, mesh.GetNodeAt(i)->GetPoint()[0], i);
////            mpNewTestWriter->PutVariable(new_p_var_id, cellPressureAns[i], i);
////            mpNewTestWriter->PutVariable(new_x_var_id, mesh.GetNodeAt(100)->GetPoint()[0], i);
////            mpNewTestWriter->PutVariable(new_c_var_id, concentrationAns[i], i); 
////        }
////        mpNewTestWriter->Close();
////   
////        //read in good data file and compare line by line
////        std::ifstream testfile("data/CancerQ5.dat",std::ios::in);
////        std::ifstream goodfile("data/CancerQ5Good.dat",std::ios::in);
////        std::string teststring;
////        std::string goodstring;
////        while(getline(testfile, teststring))
////        {
////              getline(goodfile,goodstring);
////              TS_ASSERT_EQUALS(teststring,goodstring);
////        }
////        testfile.close();
////        goodfile.close();
//        
//    }
//

};

#endif //_TESTPRACTICALONE_HPP_
