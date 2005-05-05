#ifndef _TESTPRACTICALTWO_HPP_
#define _TESTPRACTICALTWO_HPP_

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

#include "ConcentrationPdePrac2.hpp"
#include "CellPressurePdePrac2.hpp"
#include "BoundaryOdeSystem.hpp"

#include "SimpleLinearEllipticAssembler.hpp"
#include "SimpleNonlinearEllipticAssembler.hpp"

#include "EulerIvpOdeSolver.hpp"

  
class TestPracticalTwo : public CxxTest::TestSuite 
{
public:


    void testPrac2Question1(void)
    {
        int FakeArgc=0;
        char *FakeArgv0="testrunner";
        char **FakeArgv=&FakeArgv0;
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
        
        const int SPACE_DIM = 1;
        const int ELEMENT_DIM = 1;
        double alpha = 0.5;
        double rho = 1.0;
        int endTimeIndex = 1;
        double time = 0.0;
        double timestep = 0.01;
        double temp;
        double X_new = 1.0;
        double tol = 0.0001;
                
         // Create mesh from mesh reader
        
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
        ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        int numberOfNodes = mesh.GetNumNodes();
        int numberOfElements = mesh.GetNumElements();
                   
        // Instantiate PDE object for oxygen concentration
        ConcentrationPdePrac2<SPACE_DIM> concentrationPde;
        concentrationPde.mNumberOfNodes=numberOfNodes;
        
        // Instantiate PDE object for cell Pressure
        CellPressurePdePrac2<SPACE_DIM> cellPressurePde;  
        cellPressurePde.mRho = rho;
        
        // Boundary conditions for oxygen concentration
        // C'(0)=0, C(1)=1
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM> concentrationBcc;
        ConstBoundaryCondition<SPACE_DIM>* pConcentrationBoundaryCondition1 = new ConstBoundaryCondition<SPACE_DIM>(0.0);
        
        ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator concentrationIter = mesh.GetFirstBoundaryElement();
        concentrationBcc.AddNeumannBoundaryCondition(*concentrationIter,pConcentrationBoundaryCondition1);
        
        ConstBoundaryCondition<SPACE_DIM>* pConcentrationBoundaryCondition2 = new ConstBoundaryCondition<SPACE_DIM>(1.0);
        concentrationBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pConcentrationBoundaryCondition2);
        
        // Boundary conditions for cell pressure
        // Pc'(0)=0, Pc(1)=0
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM> cellPressureBcc;
        ConstBoundaryCondition<SPACE_DIM>* pCellPressureBoundaryCondition1 = new ConstBoundaryCondition<SPACE_DIM>(0.0);
        
        ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator cellPressureIter = mesh.GetFirstBoundaryElement();
        cellPressureBcc.AddNeumannBoundaryCondition(*cellPressureIter,pCellPressureBoundaryCondition1);
        
        // RHS cellPressure BC taken into time loop as it is now time dependent
        
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
        SimpleNonlinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM> concentrationAssembler;
        
        // Assembler for cell pressure
        SimpleLinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM> cellPressureAssembler;
        
        
        // Set up solution guess for residuals
        int length=numberOfNodes;
        
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
        
        GaussianQuadratureRule<ELEMENT_DIM> quadRule(2);
        LinearBasisFunction<ELEMENT_DIM> basis_func;
        
        Vec concentrationAnswer;
        Vec concentrationResidual;
        VecDuplicate(initialGuess,&concentrationResidual);
        VecDuplicate(initialGuess,&concentrationAnswer);
        
        
        // Set up things for use in time loop...
        double scaleFactor = 0.0;
        double pressureAtBoundary = 0.0;
        double *cellPressureAns;
        double *concentrationAns;
        double cellVelocity[numberOfNodes];
        double xNecrotic = 0.0;
        double x_i = 0.0;
        double x_iMinus1 = 0.0;
        int i=0;
        double x_alpha = 0.0;
        int cellPressureIerr;
        double fluidPressureAns[numberOfNodes];
        int concentrationIerr;
        const double phi0 = 0.5;
        double fluidVelocity[numberOfNodes];
        // Set up the phi vector...
        std::vector<double> phi;
        for(int j = 0 ; j < numberOfNodes ; j++ )
        {
            phi.push_back(phi0);
        } 
        // Start of grand time loop
        for(int timeIndex = 0 ; timeIndex < endTimeIndex ; timeIndex++)
        {
            // Define pressureAtBoundary, p(t), where p(t)=Pc + Pe.
            pressureAtBoundary = 1.0;
            concentrationPde.mPhi = phi;
            concentrationPde.mXnecrotic = xNecrotic;
            concentrationPde.mX = X_new;
            try {
                concentrationAnswer=concentrationAssembler.AssembleSystem(&mesh, &concentrationPde, &concentrationBcc, &concentrationSolver, &basis_func, &quadRule, initialGuess, true);
            } catch (Exception e) {
                TS_TRACE(e.getMessage());
            }
                                          
            // Check result
            
            concentrationIerr = VecGetArray(concentrationAnswer, &concentrationAns);
            
            // The Solver fails if the initial guess is updated!!
            
            // Calculate the value of x at which C = alpha , this boundary determines where the cells begin to die of oxygen deprivation
            
            i = 0;
            x_alpha = 0.0;
            while(concentrationAns[i] < alpha && i<=numberOfNodes)
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
            
            // Just this BC in the time loop, because it is time dependent
            ConstBoundaryCondition<SPACE_DIM>* pCellPressureBoundaryCondition2 = new ConstBoundaryCondition<SPACE_DIM>(pressureAtBoundary);
            cellPressureBcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pCellPressureBoundaryCondition2);
            
            // Solve the linear pde for cell Pressure                
            Vec cellPressureAnswer = cellPressureAssembler.AssembleSystem(mesh, &cellPressurePde, cellPressureBcc, &cellPressureSolver);
            
//            // Test to see if result lies within bounds calculated from Matlab
            
            cellPressureIerr = VecGetArray(cellPressureAnswer, &cellPressureAns);
          //Solution should lie between 0.06 and -0.07
            for (int j=0; j < numberOfNodes; j++)
            {
//                std::cout << "Cell Pressure " << j << " = " << cellPressureAns[j] << "\n";
//                TS_ASSERT_DELTA(cellPressureAns[j], 0.07, 0.14);
//                TS_ASSERT_DELTA(cellPressureAns[j], -0.07, 0.14); 
            }
            VecRestoreArray(cellPressureAnswer, &cellPressureAns);
            
            
            // Solve for extracellular fluid pressure
            
            
            for (int j=0; j < numberOfNodes; j++)
            {
                fluidPressureAns[j] = pressureAtBoundary - cellPressureAns[j];
            }
            

            
//            // Solution should lie between -0.06 and 0.07 ,bounds obtained from matlab
//            for (int i=0; i < numberOfNodes; i++)
//            {
//                TS_ASSERT_DELTA(fluidPressureAns[i], -0.07, 0.14);
//                TS_ASSERT_DELTA(fluidPressureAns[i], 0.07, 0.14); 
//            }
            
            // Solve for cell velocity , U_c partial derivative of cell pressure w.r.t x
            
            
            
                       
            for (int j=0; j < numberOfNodes; j++)
            {
                if (j == 0)
                {
                    // Use forward difference
                    cellVelocity[j] = - (cellPressureAns[j+1] - cellPressureAns[j])*numberOfElements;
                }
                else if(j == numberOfElements)
                {
                    // Use backward difference
                    cellVelocity[j] = - (cellPressureAns[j] - cellPressureAns[j-1])*numberOfElements;
                }
                else
                {
                    // Use central difference
                    cellVelocity[j] = - (cellPressureAns[j+1] - cellPressureAns[j-1])*numberOfElements/2.0;
                }
                
            }
            
//            //Test cell pressure lies within bounds from matlab
//            for (int j=0; j < numberOfNodes; j++)
//            {
//                TS_ASSERT_DELTA(cellPressureAns[j], -0.4, 0.70001);
//                TS_ASSERT_DELTA(cellPressureAns[j], 0.3, 0.70001); 
//            }
                  
            // Solve for U_e

            for (int j=0; j < numberOfNodes; j++)
            {
                fluidVelocity[j] = (-phi0 / (1.0 - phi0))*cellVelocity[j];
            }
//            for (int i=0; i < numberOfNodes; i++)
//            {
//                TS_ASSERT_DELTA(cellPressureAns[i], -0.3, 0.70001);
//                TS_ASSERT_DELTA(cellPressureAns[i], 0.4, 0.70001); 
//            }

            std::cout << "Cell Velocity at X = " << cellVelocity[numberOfElements] << "\n";                         
            pBoundarySystem->mCellVelocityAtBoundary = cellVelocity[numberOfElements];
            
            solutions = myEulerSolver->Solve(pBoundarySystem, time, time+timestep, timestep, boundaryInit);  
            // Update xNecrotic for next time step
            
            // Find where cellPressure[i] = fluidPressure[i] ... if it doesn't happen then Xnecrotic=0
            i = 0;
            xNecrotic = 0.0;
            while(cellPressureAns[i] - fluidPressureAns[i] < tol)
            {
                i++;
            }
            
            if (i > 0) 
            {
                x_i = mesh.GetNodeAt(i)->GetPoint()[0];
                x_iMinus1 = mesh.GetNodeAt(i-1)->GetPoint()[0];
                xNecrotic = x_iMinus1 + (x_i - x_iMinus1)*(pressureAtBoundary/2.0 - cellPressureAns[i-1])/(cellPressureAns[i]-cellPressureAns[i-1]);
            }
            cellPressurePde.mXnecrotic = xNecrotic;
            
            for(int j = 0 ; j < i ; j++)
            {
                cellPressureAns[j]=pressureAtBoundary/2.0;
                fluidPressureAns[j]=pressureAtBoundary/2.0;
                // in here reduce all phi s by a factor of exp-rho timestep.
                phi[j]=phi[j]*exp(-rho*timestep);
                cellVelocity[j]=0.0;
                fluidVelocity[j]=0.0;
            }
            
                     
            // Rescale Mesh
            // Get X_new, rescale with X_old
            
            X_new = solutions.mSolutions[1][0];
            std::cout << "X_new = "<< X_new << std::endl;
            scaleFactor = solutions.mSolutions[1][0] / solutions.mSolutions[0][0];
            //std::cout << "Scale Factor = " << scaleFactor << std::endl;
            
            for(int k=0; k<numberOfNodes; k++)
            {
                for(int j=0; j<SPACE_DIM ; j++)
                {
                    temp = scaleFactor*mesh.GetNodeAt(k)->GetPoint()[j];
                    //std::cout << "new location = " << temp << std::endl;
                    mesh.GetNodeAt(k)->GetPoint().SetCoordinate(j,temp);
                }
            }
            //std::cout << "New RHS node location = " << mesh.GetNodeAt(numberOfNodes-1)->GetPoint()[0] << std::endl;           
            time = time + timestep;
        // End of grand time loop
        }    
        
        // set up file output
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","CancerPrac2Q1");
        mpNewTestWriter->DefineFixedDimension("Space","dimensionless", numberOfNodes);
        int new_c_var_id = mpNewTestWriter->DefineVariable("Concentration","mM");
        int new_space_var_id = mpNewTestWriter->DefineVariable("Space","dimensionless");
        int new_x_var_id = mpNewTestWriter->DefineVariable("X","dimensionless");
        int new_p_var_id = mpNewTestWriter->DefineVariable("Pressure","Pa");
        int cellV_var_id = mpNewTestWriter->DefineVariable("Cell Velocity","dimensionless");
        mpNewTestWriter->EndDefineMode();
        
        
        for (int i = 0; i < numberOfNodes; i++) 
        {
                 
            mpNewTestWriter->PutVariable(new_space_var_id, mesh.GetNodeAt(i)->GetPoint()[0], i);
            mpNewTestWriter->PutVariable(new_p_var_id, cellPressureAns[i], i);
            mpNewTestWriter->PutVariable(new_x_var_id, mesh.GetNodeAt(100)->GetPoint()[0], i);
            mpNewTestWriter->PutVariable(new_c_var_id, concentrationAns[i], i);
            mpNewTestWriter->PutVariable(cellV_var_id, cellVelocity[i], i);
        }
        mpNewTestWriter->Close();
//   
//        //read in good data file and compare line by line
//        std::ifstream testfile("data/CancerQ5.dat",std::ios::in);
//        std::ifstream goodfile("data/CancerQ5Good.dat",std::ios::in);
//        std::string teststring;
//        std::string goodstring;
//        while(getline(testfile, teststring))
//        {
//              getline(goodfile,goodstring);
//              TS_ASSERT_EQUALS(teststring,goodstring);
//        }
//        testfile.close();
//        goodfile.close();
        
    }


};

#endif //_TESTPRACTICALTWO_HPP_
