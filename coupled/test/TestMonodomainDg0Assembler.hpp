#ifndef _TESTMONODOMAINDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINDG0ASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"  
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "PropagationPropertiesCalculator.hpp"

#include "MonodomainPde.hpp"
#include "MockEulerIvpOdeSolver.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractMonodomainProblemStimulus.hpp"

#include "FakeLr91Pde.hpp"

class PointStimulus1D: public AbstractMonodomainProblemStimulus<1>
{
public:
    virtual void Apply(MonodomainPde<1> *pPde, 
                       ConformingTetrahedralMesh<1,1> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);

        for (int i = 0; i < pMesh->GetNumNodes(); i++)
        {
            if (pMesh->GetNodeAt(i)->GetPoint()[0] == 0.0)
            {
                pPde->SetStimulusFunctionAtNode(i, &stimulus);
            }
        }
    }
};

class EdgeStimulus2D: public AbstractMonodomainProblemStimulus<2>
{
    virtual void Apply(MonodomainPde<2> *pPde,
                       ConformingTetrahedralMesh<2,2> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);
        
        for (int i = 0; i < pMesh->GetNumNodes(); i++)
        {
            if (pMesh->GetNodeAt(i)->GetPoint()[0] == 0.0)
            {
                pPde->SetStimulusFunctionAtNode(i, &stimulus);
            }
        }
    }
};

class PointStimulus2D: public AbstractMonodomainProblemStimulus<2>
{
private:
    int mNode;
    
public:
    PointStimulus2D(const int node = 0)
    {
        mNode = node;
    }
    
    virtual void Apply(MonodomainPde<2> *pPde,
                       ConformingTetrahedralMesh<2,2> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);

        pPde->SetStimulusFunctionAtNode(mNode, &stimulus);
    }
};

class FaceStimulus3D: public AbstractMonodomainProblemStimulus<3>
{
    virtual void Apply(MonodomainPde<3> *pPde,
                       ConformingTetrahedralMesh<3,3> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);

        for (int i = 0; i < pMesh->GetNumNodes(); i++)
        {
            if (pMesh->GetNodeAt(i)->GetPoint()[0] == 0.0)
            {
                pPde->SetStimulusFunctionAtNode(i, &stimulus);
            }
        }
    }
};

class PointStimulus3D: public AbstractMonodomainProblemStimulus<3>
{
    virtual void Apply(MonodomainPde<3> *pPde,
                       ConformingTetrahedralMesh<3,3> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);

        pPde->SetStimulusFunctionAtNode(0, &stimulus);
    }
};

class TestMonodomainDg0Assembler : public CxxTest::TestSuite 
{   
private:
	/**
	 * Refactor code to set up a PETSc vector holding the initial condition.
	 */
	Vec CreateInitialConditionVec(int size)
	{
    	Vec initial_condition;
    	VecCreate(PETSC_COMM_WORLD, &initial_condition);
    	VecSetSizes(initial_condition, PETSC_DECIDE, size);
    	VecSetFromOptions(initial_condition);
    	return initial_condition;
	}
    
public:
	void TestMonodomainDg01D()
	{
        PointStimulus1D point_stimulus_1D;
        MonodomainProblem<1> monodomainProblem("mesh/test/data/1D_0_to_1_100_elements",
                                               30, 
                                               "testoutput/MonoDg01d",
                                               "NewMonodomainLR91_1d",
                                               &point_stimulus_1D);

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }   
            }
               
            if (global_index==1)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], 9.4435, 0.001);
            }
            if (global_index==3)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], 9.4145, 0.001);
            }
            if (global_index==5)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], 9.3630, 0.001);
            }
            if (global_index==7)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], 9.2948, 0.001);
            }
            if (global_index==9)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], 9.2135, 0.001);
            }
            if (global_index==75)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], 5.4796, 0.001);
            }
        }





        // Display the voltage at the right hand node of the mesh.
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
           for (int i = 0; i < monodomainProblem.mMesh.GetNumNodes(); i++)
            {
                if (monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[0] == 1.0)
                {
//std::cout << "Voltage: " << currentVoltageArray[i] << std::endl;
                }
            }
        }





        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);

    }
    
    void TestMonodomainDg02DWithEdgeStimulus( void )
    {   
        EdgeStimulus2D edge_stimulus_2D;
        
        MonodomainProblem<2> monodomainProblem("mesh/test/data/square_128_elements",
                                               0.01,//0.1,
// 0.01 is just so that we run the model for 1 time step only. will change to match prev test
                                               "testoutput/MonoDg02dWithEdgeStimulus",
                                               "NewMonodomainLR91_2dWithEdgeStimulus",
                                               &edge_stimulus_2D);

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
            }
        }

        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            /*
             * Test the top right node against the right one in the 1D case, 
             * comparing voltage and then test all the nodes on the right hand 
             * side of the square against the top right one, comparing voltage
             */

            bool needInitialisation = true;
            double rhsVoltage;

//std::cout << "--------------------------------------------------" << std::endl;
//std::cout << "LEFT HAND SIDE" << std::endl;
//
//            for (int i = 0; i < monodomainProblem.mMesh.GetNumNodes(); i++)
//            {
//                if (monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[0] == 0.0)
//                {
//                    // x = 0 is where the stimulus has been applied
//                    // x = 1 is the other end of the mesh and where we want to 
//                    //       to test the value of the nodes
//                    
//                    if (needInitialisation)
//                    {
//                        rhsVoltage = currentVoltageArray[i];
//                        needInitialisation = false;
//
////                        TS_ASSERT_DELTA(rhsVoltage, -84.5001, 0.001);   // No diffusion (tEnd = 0.01 ms)
////                        TS_ASSERT_DELTA(rhsVoltage, -84.5015, 0.001);   // No diffusion (tEnd = 0.1 ms)
//                    }
//                    else
//                    {
//                        TS_ASSERT_DELTA(currentVoltageArray[i], rhsVoltage, 0.001);
//                        std::cout << "y=" << monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
//                    }
//                }
//            }

//std::cout << "--------------------------------------------------" << std::endl;
//std::cout << "RIGHT HAND SIDE" << std::endl;

            needInitialisation = true;

            for (int i = 0; i < monodomainProblem.mMesh.GetNumNodes(); i++)
            {
                if (monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[0] == 1.0)
                {
                    // x = 0 is where the stimulus has been applied
                    // x = 1 is the other end of the mesh and where we want to 
                    //       to test the value of the nodes
                    
                    if (needInitialisation)
                    {
                        rhsVoltage = currentVoltageArray[i];
                        needInitialisation = false;

//                        TS_ASSERT_DELTA(rhsVoltage, -84.5001, 0.001);   // No diffusion (tEnd = 0.01 ms)
//                        TS_ASSERT_DELTA(rhsVoltage, -84.5015, 0.001);   // No diffusion (tEnd = 0.1 ms)
                    }
                    else
                    {
                        TS_ASSERT_DELTA(currentVoltageArray[i], rhsVoltage, 0.001);
                       // std::cout << "y=" << monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
                    }
                }
            }
        }
        
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   

    void TestMonodomainDg02DWithPointStimulus( void )
    {   
        PointStimulus2D point_stimulus_2D(4); // Central node
        
        MonodomainProblem<2> monodomainProblem("mesh/test/data/square_128_elements",
                                               0.1, 
                                               "testoutput/MonoDg02dWithPointStimulus",
                                               "NewMonodomainLR91_2dWithPointStimulus",
                                               &point_stimulus_2D);

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
            }
        }
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            /*
             * Test that corners are equal, and centres of sides.
             */
            TS_ASSERT_DELTA(currentVoltageArray[0], currentVoltageArray[1], 0.0001);
            TS_ASSERT_DELTA(currentVoltageArray[0], currentVoltageArray[2], 0.0001);
            TS_ASSERT_DELTA(currentVoltageArray[0], currentVoltageArray[3], 0.0001);
            
            TS_ASSERT_DELTA(currentVoltageArray[5], currentVoltageArray[6], 0.0001);
            TS_ASSERT_DELTA(currentVoltageArray[5], currentVoltageArray[7], 0.0001);
            TS_ASSERT_DELTA(currentVoltageArray[5], currentVoltageArray[8], 0.0001);
        }
        
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   

    void TestMonodomainDg03DWithFaceStimulus( void )
    {
        FaceStimulus3D face_stimulus_3D;
        
        MonodomainProblem<3> monodomainProblem("mesh/test/data/slab_395_elements",
                                               0.25, 
                                               "testoutput/MonoDg03dWithFaceStimulus",
                                               "NewMonodomainLR91_3dWithFaceStimulus",
                                               &face_stimulus_3D);
        ///todo This test now fails (not sensible ODE values) for mEndTime >= 0.26 

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
            }
        }
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   

    void TestMonodomainDg03DWithPointStimulus( void )
    {
        PointStimulus3D point_stimulus_3D;
        
        MonodomainProblem<3> monodomainProblem("mesh/test/data/slab_395_elements",
                                               0.25, 
                                               "testoutput/MonoDg03dWithPointStimulus",
                                               "NewMonodomainLR91_3dWithPointStimulus",
                                               &point_stimulus_3D);
        ///todo This test now fails (not sensible ODE values) for mEndTime >= 0.26 

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
            }
        }
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }
    
    
    
    void TestFakePdeWithSimpleAssembler2D( void )
    {       
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        FakeLr91Pde<2> pde;         
    
        // Boundary conditions - zero neumann on boundary;
        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
//        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
//        ConstBoundaryCondition<2>* zero_neumann_bc = new ConstBoundaryCondition<2>(0.0);        
//        while(iter < mesh.GetLastBoundaryElement())
//        {
//            bcc.AddNeumannBoundaryCondition(*iter, zero_neumann_bc);
//            iter++;
//        }
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetFirstBoundaryNode();
        ConstBoundaryCondition<2>* zero_dirichlet_bc = new ConstBoundaryCondition<2>(-84.5);        
        while(iter < mesh.GetLastBoundaryNode())
        {
            bcc.AddDirichletBoundaryCondition(*iter, zero_dirichlet_bc);
            iter++;
        }           
        // Linear solver
        SimpleLinearSolver linearSolver;
    
        // Assembler
        SimpleDg0ParabolicAssembler<2,2> fullSolver;
        
        // initial condition;   
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, mesh.GetNumNodes() );
        VecSetFromOptions(initial_condition);
  
        double* initial_condition_array;
        VecGetArray(initial_condition, &initial_condition_array); 
        
        int mLo, mHi;
        VecGetOwnershipRange(initial_condition, &mLo, &mHi);
        
        // Set a constant initial voltage throughout the mMesh
        for(int global_index=mLo; global_index<mHi; global_index++)
        {
            initial_condition_array[global_index-mLo] = -84.5;
        }
        VecRestoreArray(initial_condition, &initial_condition_array);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);
        
        double t_end = 0.01;
        fullSolver.SetTimes(0, t_end, 0.01);
        fullSolver.SetInitialCondition(initial_condition);

        Vec result = fullSolver.Solve(mesh, &pde, bcc, &linearSolver);
        
        double* currentVoltageArray;
        int ierr = VecGetArray(result, &currentVoltageArray); 
 
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(currentVoltageArray[i], -84.5, 0.0001);
//                std::cout << "x=" << mesh.GetNodeAt(i)->GetPoint()[0] << std::endl;
//                std::cout << "y=" << mesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
//                
            }
        }
        
        VecRestoreArray(result, &currentVoltageArray);
        
        VecDestroy(initial_condition);
        VecDestroy(result);
    }
    
    void TestFakePdeWithSimpleAssembler1D( void )
    {       
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        FakeLr91Pde<1> pde;
    
        // Boundary conditions - zero neumann on boundary;
        BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
//        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
//        ConstBoundaryCondition<2>* zero_neumann_bc = new ConstBoundaryCondition<2>(0.0);        
//        while(iter < mesh.GetLastBoundaryElement())
//        {
//            bcc.AddNeumannBoundaryCondition(*iter, zero_neumann_bc);
//            iter++;
//        }
        ConformingTetrahedralMesh<1,1>::BoundaryNodeIterator iter = mesh.GetFirstBoundaryNode();
        ConstBoundaryCondition<1>* zero_dirichlet_bc = new ConstBoundaryCondition<1>(-84.5);
        while(iter < mesh.GetLastBoundaryNode())
        {
            bcc.AddDirichletBoundaryCondition(*iter, zero_dirichlet_bc);
            iter++;
        }           
        // Linear solver
        SimpleLinearSolver linearSolver;
    
        // Assembler
        SimpleDg0ParabolicAssembler<1,1> fullSolver;
        
        // initial condition;   
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, mesh.GetNumNodes() );
        VecSetFromOptions(initial_condition);
  
        double* initial_condition_array;
        VecGetArray(initial_condition, &initial_condition_array); 
        
        int mLo, mHi;
        VecGetOwnershipRange(initial_condition, &mLo, &mHi);
        
        // Set a constant initial voltage throughout the mMesh
        for(int global_index=mLo; global_index<mHi; global_index++)
        {
            initial_condition_array[global_index-mLo] = -84.5;
        }
        VecRestoreArray(initial_condition, &initial_condition_array);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);
        
        double t_end = 1;
        fullSolver.SetTimes(0, t_end, 0.01);
        fullSolver.SetInitialCondition(initial_condition);

        Vec result = fullSolver.Solve(mesh, &pde, bcc, &linearSolver);
        
        double* currentVoltageArray;
        int ierr = VecGetArray(result, &currentVoltageArray); 
 
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(currentVoltageArray[i], -84.5, 0.0001);
  //              std::cout << "x=" << mesh.GetNodeAt(i)->GetPoint()[0] << std::endl;
                
            }
        }
        
        VecRestoreArray(result, &currentVoltageArray);
        
        VecDestroy(initial_condition);
        VecDestroy(result);
    }
};
#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
