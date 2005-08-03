#ifndef _TESTPARALLELMONODOMAINDG0ASSEMBLER_HPP_
#define _TESTPARALLELMONODOMAINDG0ASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <iostream>
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"  
#include "TrianglesMeshReader.hpp"
#include "ParallelMonodomainPde.hpp"
#include "ParallelMonodomainDg0Assembler.hpp"
#include "ColumnDataWriter.hpp"
#include <cmath>
#include "PetscSetupAndFinalize.hpp"

 
class TestMonodomainDg0Assembler : public CxxTest::TestSuite 
{   
public:
	void testMonodomainDg01D()
	{  
        
		
		double tStart = 0; 
		double tFinal = 0.1;
        
		// use big time step (the pde timestep) is the same as the small time step (the ode timestep)
		double tBigStep = 0.01; 
		double tSmallStep  = 0.01;
        
		// Create mesh from mesh reader 
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
        
		// Instantiate PDE object
		AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
		ParallelMonodomainPde<1> monodomain_pde(mesh.GetNumNodes(), pMySolver, tStart, tBigStep, tSmallStep);
        
		// sets Luo Rudy system with initial conditions passed on
		double voltage = -9999; // This voltage will be ignored
		double m = 0.0017;
		double h = 0.9833;
		double j = 0.9895;  
		double d = 0.003;
		double f = 1;
		double x = 0.0056;
		double caI = 0.0002;
		double magnitudeOfStimulus = -80.0;  
		double durationOfStimulus  = 0.5 ;
		          
		// bad 
		std::vector<double> initialConditions;
		initialConditions.push_back(h);
		initialConditions.push_back(j);
		initialConditions.push_back(m);
		initialConditions.push_back(caI);
		initialConditions.push_back(voltage);
		initialConditions.push_back(d);
		initialConditions.push_back(f);
		initialConditions.push_back(x);

		// set this as the initial condition of the gating vars at each node in the mesh        
		monodomain_pde.SetUniversalInitialConditions(initialConditions);
        
		//monodomain_pde.SetInitialConditions(mesh, 
		//need to write mesh.GetInitialConditionsNodes
        
		// add initial stim to node 0 only
		AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus);
		monodomain_pde.SetStimulusFunctionAtNode(0, pStimulus);
                
        
		// Boundary conditions: zero neumann on entire boundary
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		ConstBoundaryCondition<1>* pNeumannBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
		ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
		bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition1);

		ConstBoundaryCondition<1>* pNeumannBoundaryCondition2 = new ConstBoundaryCondition<1>(0.0);
		iter = mesh.GetLastBoundaryElement();
		iter--;
		bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition2);
        
		// Linear solver
		SimpleLinearSolver linearSolver;
    
		// Assembler
		ParallelMonodomainDg0Assembler<1,1> monodomainAssembler;
        
		// initial condition;   
		Vec currentVoltage;
		VecCreate(PETSC_COMM_WORLD, &currentVoltage);
		VecSetSizes(currentVoltage, PETSC_DECIDE, mesh.GetNumNodes() );
		//VecSetType(initialCondition, VECSEQ);
		VecSetFromOptions(currentVoltage);
  
		double* currentVoltageArray;
		int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
		// initial voltage condition of a constant everywhere on the mesh
		int lo, hi;
		VecGetOwnershipRange(currentVoltage,&lo,&hi);
		for(int i=0; i<hi-lo; i++)
		{
			currentVoltageArray[i] = -84.5;
		}
		VecRestoreArray(currentVoltage, &currentVoltageArray);      
		VecAssemblyBegin(currentVoltage);
		VecAssemblyEnd(currentVoltage);

		// Hackish fix to save typing...
		Vec initial_condition = currentVoltage;
	          
		/*
		 * Write data to a file NewMonodomainLR91_1d_xx.dat, 'xx' refers to nth time step
		 *  using ColumnDataWriter 
		 */                                                                            
           
        
        // uncomment all column writer related lines to write data (and the line further below)         
        //ColumnDataWriter *mpTestWriter;
        //mpTestWriter = new ColumnDataWriter("data","NewMonodomainLR91_1d");
       
      //  int time_var_id = 0;
       // int voltage_var_id = 0;

       // mpTestWriter->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
       // mpTestWriter->DefineUnlimitedDimension("Time","msecs");

      //  time_var_id = mpTestWriter->DefineVariable("Time","msecs");
     //   voltage_var_id = mpTestWriter->DefineVariable("V","mV");
     //   mpTestWriter->EndDefineMode();
           
		double tCurrent = tStart;        
		while( tCurrent < tFinal )
        { 
            // std::cout << "t = " << tCurrent << "\n" << std::flush;
	
            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( initial_condition );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
    
    		// New initial condition is current solution
    		VecDestroy(initial_condition);
    		initial_condition = currentVoltage;
	        
            // Writing data out to the file NewMonodomainLR91_1d.dat
         
        //    int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
              
         //   mpTestWriter->PutVariable(time_var_id, tCurrent); 
          //  for(int j=0; j<mesh.GetNumNodes(); j++) 
          //  {
               // uncomment the line below to write data (and the line further below)
               // mpTestWriter->PutVariable(voltage_var_id, currentVoltageArray[j], j);    
         //   }
  
        //    VecRestoreArray(currentVoltage, &currentVoltageArray); 
            // uncomment the line below to write data
            // mpTestWriter->AdvanceAlongUnlimitedDimension();
     
            monodomain_pde.ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
        }

        // close the file that stores voltage values
        //mpTestWriter->Close();

        // test whether voltages and gating variables are in correct ranges
        ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        //int lo,hi;
        //VecGetOwnershipRange(currentVoltage,&lo,&hi);
        for(int globalIndex=lo; globalIndex<hi; globalIndex++)
        {
           
            	// assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[globalIndex-lo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[globalIndex-lo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_pde.GetOdeVarsAtNode(globalIndex);           
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
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(currentVoltage);
        VecAssemblyEnd(currentVoltage);
        
        delete pMySolver;
        delete pStimulus;
        VecDestroy(currentVoltage);
    }
    
 
    
    void testMonodomainDg02D( void )
    {   
        double tStart = 0; 
        double tFinal = 0.1;
        
        // use big time step (the pde timestep) is the same as the small time step (the ode timestep)
        double tBigStep = 0.01; 
        double tSmallStep  = 0.01;
        
        // read mesh on [0,1]x[0,1]
        TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        ParallelMonodomainPde<2> monodomain_pde(mesh.GetNumNodes(), pMySolver, tStart, tBigStep, tSmallStep);
                
        
        // sets Luo Rudy system with initial conditions passed on
        double voltage = -9999; // This voltage will be ignored
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;
        double d = 0.003;
        double f = 1;
        double x = 0.0056;
        double caI = 0.0002;
        double magnitudeOfStimulus = -80.0;  
        double durationOfStimulus  = 0.5 ;
                  
        // bad 
        std::vector<double> initialConditions;
        initialConditions.push_back(h);
        initialConditions.push_back(j);
        initialConditions.push_back(m);
        initialConditions.push_back(caI);
        initialConditions.push_back(voltage);
        initialConditions.push_back(d);
        initialConditions.push_back(f);
        initialConditions.push_back(x);

        // set this as the initial condition of the gating vars at each node in the mesh        
        monodomain_pde.SetUniversalInitialConditions(initialConditions);
                
        // add initial stim to node 0 only
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus);
        monodomain_pde.SetStimulusFunctionAtNode(0, pStimulus);
                         
        // Boundary conditions: zero neumann on boundary
        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetFirstBoundaryElement();
        ConstBoundaryCondition<2>* pNeumannBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        
        while(surf_iter < mesh.GetLastBoundaryElement())
        {
            bcc.AddNeumannBoundaryCondition(*surf_iter, pNeumannBoundaryCondition);
            surf_iter++;
        }
        
        
        // Linear solver
        SimpleLinearSolver linearSolver;
    
        // Assembler
        ParallelMonodomainDg0Assembler<2,2> monodomainAssembler;
        
        // initial condition;   
        Vec currentVoltage;
        VecCreate(PETSC_COMM_WORLD, &currentVoltage);
        VecSetSizes(currentVoltage, PETSC_DECIDE, mesh.GetNumNodes() );
        //VecSetType(initialCondition, VECSEQ);
        VecSetFromOptions(currentVoltage);
  
        double* currentVoltageArray;
        int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        // initial voltage condition of a constant everywhere on the mesh
        int lo, hi;
        VecGetOwnershipRange(currentVoltage,&lo,&hi);
        
        for(int i=0; i<hi-lo; i++)
        {
            currentVoltageArray[i] = -84.5;
        }
     
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(currentVoltage);
        VecAssemblyEnd(currentVoltage);

		// Hackish fix to save typing...
		Vec initial_condition = currentVoltage;
              
        /*
        * Write data to a file NewMonodomainLR91_2d_xx.dat, 'xx' refers to nth time step
        *  using ColumnDataWriter 
        */                                                           
                 
     
        
        // uncomment all column writer related lines to write data (and the line further below)         
        //ColumnDataWriter *mpTestWriter;
        //mpTestWriter = new ColumnDataWriter("data","NewMonodomainLR91_2d");
       
      //  int time_var_id = 0;
       // int voltage_var_id = 0;

       // mpTestWriter->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
       // mpTestWriter->DefineUnlimitedDimension("Time","msecs");

      //  time_var_id = mpTestWriter->DefineVariable("Time","msecs");
     //   voltage_var_id = mpTestWriter->DefineVariable("V","mV");
     //   mpTestWriter->EndDefineMode();
           
        double tCurrent = tStart;        
        while( tCurrent < tFinal )
        {
            // std::cout << "t = " << tCurrent << "\n" << std::flush;

            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( initial_condition );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            
            // New initial condition is current solution
    		VecDestroy(initial_condition);
    		initial_condition = currentVoltage;
            
            // Writing data out to the file NewMonodomainLR91_2d.dat
         
        //    int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
              
         //   mpTestWriter->PutVariable(time_var_id, tCurrent); 
          //  for(int j=0; j<mesh.GetNumNodes(); j++) 
          //  {
               // uncomment the line below to write data (and the line further below)
               // mpTestWriter->PutVariable(voltage_var_id, currentVoltageArray[j], j);    
         //   }
  
        //    VecRestoreArray(currentVoltage, &currentVoltageArray); 
            // uncomment the line below to write data
            // mpTestWriter->AdvanceAlongUnlimitedDimension();
     
            monodomain_pde.ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
        }

        // close the file that stores voltage values
        //mpTestWriter->Close();
    
            
        // this saves to one file all voltages at one chosen node, rather than creating a series of files      
        /*         
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("data","NewMonodomainLR91_2d");
        mpNewTestWriter->DefineUnlimitedDimension("Time","msecs");
        int new_time_var_id = mpNewTestWriter->DefineVariable("Time","msecs");
        int new_v_var_id = mpNewTestWriter->DefineVariable("V","mV");   
        mpNewTestWriter->EndDefineMode();
           
        double tCurrent = tStart;        
        while( tCurrent < tFinal )
        {
            // std::cout << "t = " << tCurrent << "\n" << std::flush;
            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( currentVoltage );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            
            // Writing data out to the file NewMonodomainLR91.dat
         
            int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
              
            // output to file
            mpNewTestWriter->PutVariable(new_time_var_id, tCurrent);
            mpNewTestWriter->PutVariable(new_v_var_id, currentVoltageArray[53]);
            VecRestoreArray(currentVoltage, &currentVoltageArray); 

            mpNewTestWriter->AdvanceAlongUnlimitedDimension();
              
            monodomain_pde.ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
        }
        
        mpNewTestWriter->Close();
        */



        // test whether voltages and gating variables are in correct ranges
        ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        for(int globalIndex=lo; globalIndex<hi; globalIndex++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[globalIndex-lo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[globalIndex-lo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_pde.GetOdeVarsAtNode(globalIndex);           
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
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(currentVoltage);
        VecAssemblyEnd(currentVoltage);
 
		delete pMySolver;
		delete pStimulus;
		VecDestroy(currentVoltage);
    }   


    void testMonodomainDg03D( void )
    {   
        double tStart = 0; 
        double tFinal = 0.1;
         
        
        double tBigStep = 0.01; 
        double tSmallStep  = 0.01;
        
        // read 3d mesh
        // TrianglesMeshReader mesh_reader("pdes/tests/meshdata/cylinder_with_hole_840_elements");
        TrianglesMeshReader mesh_reader("mesh/test/data/slab_138_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        
        // Instantiate PDE object
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        ParallelMonodomainPde<3> monodomain_pde(mesh.GetNumNodes(), pMySolver, tStart, tBigStep, tSmallStep);
                
        
        // sets Luo Rudy system with initial conditions passed on
        double voltage = -9999; // This voltage will be ignored
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;
        double d = 0.003;
        double f = 1;
        double x = 0.0056;
        double caI = 0.0002;
        double magnitudeOfStimulus = -80.0;  
        double durationOfStimulus  = 0.5 ;
                  
        // bad 
        std::vector<double> initialConditions;
        initialConditions.push_back(h);
        initialConditions.push_back(j);
        initialConditions.push_back(m);
        initialConditions.push_back(caI);
        initialConditions.push_back(voltage);
        initialConditions.push_back(d);
        initialConditions.push_back(f);
        initialConditions.push_back(x);
        
        // set this as the initial condition of the gating vars at each node in the mesh        
        monodomain_pde.SetUniversalInitialConditions(initialConditions);
             
        // add initial stimulus to node 0 only
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus);
        monodomain_pde.SetStimulusFunctionAtNode(0, pStimulus);
 
                 
        // Boundary conditions, zero neumann everywhere
        BoundaryConditionsContainer<3,3> bcc(1, mesh.GetNumNodes());
        ConformingTetrahedralMesh<3,3>::BoundaryElementIterator surf_iter = mesh.GetFirstBoundaryElement();
        ConstBoundaryCondition<3>* pNeumannBoundaryCondition = new ConstBoundaryCondition<3>(0.0);
        
        while(surf_iter < mesh.GetLastBoundaryElement())
        {
            bcc.AddNeumannBoundaryCondition(*surf_iter, pNeumannBoundaryCondition);
            surf_iter++;
        }
        
        // Linear solver
        SimpleLinearSolver linearSolver;
    
        // Assembler
        ParallelMonodomainDg0Assembler<3,3> monodomainAssembler;
        
        // initial condition;   
        Vec currentVoltage;
        VecCreate(PETSC_COMM_WORLD, &currentVoltage);
        VecSetSizes(currentVoltage, PETSC_DECIDE, mesh.GetNumNodes() );
        //VecSetType(initialCondition, VECSEQ);
        VecSetFromOptions(currentVoltage);
  
        double* currentVoltageArray;
        int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        // initial voltage condition of a constant everywhere on the mesh
        int lo, hi;
        VecGetOwnershipRange(currentVoltage,&lo,&hi);
        for(int i=0; i<hi-lo; i++)
        {
            currentVoltageArray[i] = -84.5;
        }
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(currentVoltage);
        VecAssemblyEnd(currentVoltage);

        // Hackish fix to save typing...
		Vec initial_condition = currentVoltage;
		
        /*
         *  Write data to a file NewMonodomainLR91_3d_xx.dat, 'xx' refers to nth time step
         *  using ColumnDataWriter 
         */                                                           
        
        
        // uncomment all column writer related lines to write data (and the line further below)         
        //ColumnDataWriter *mpTestWriter;
        //mpTestWriter = new ColumnDataWriter("data","NewMonodomainLR91_3d");
       
      //  int time_var_id = 0;
       // int voltage_var_id = 0;

       // mpTestWriter->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
       // mpTestWriter->DefineUnlimitedDimension("Time","msecs");

      //  time_var_id = mpTestWriter->DefineVariable("Time","msecs");
     //   voltage_var_id = mpTestWriter->DefineVariable("V","mV");
     //   mpTestWriter->EndDefineMode();
           
        double tCurrent = tStart;        
        while( tCurrent < tFinal )
        {
            // std::cout << "t = " << tCurrent << "\n" << std::flush;

            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( initial_condition );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            
            // New initial condition is current solution
    		VecDestroy(initial_condition);
    		initial_condition = currentVoltage;
            
            // Writing data out to the file NewMonodomainLR91_3d.dat
         
        //    int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
              
         //   mpTestWriter->PutVariable(time_var_id, tCurrent); 
          //  for(int j=0; j<mesh.GetNumNodes(); j++) 
          //  {
               // uncomment the line below to write data (and the line further below)
               // mpTestWriter->PutVariable(voltage_var_id, currentVoltageArray[j], j);    
         //   }
  
        //    VecRestoreArray(currentVoltage, &currentVoltageArray); 
            // uncomment the line below to write data
            // mpTestWriter->AdvanceAlongUnlimitedDimension();
     
            monodomain_pde.ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
        }

        // close the file that stores voltage values
        //mpTestWriter->Close();
        
        
        // test whether voltages and gating variables are in correct ranges
        ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        for(int globalIndex=lo; globalIndex<hi; globalIndex++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[globalIndex] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[globalIndex] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_pde.GetOdeVarsAtNode(globalIndex);           
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
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(currentVoltage);
        VecAssemblyEnd(currentVoltage);
        
        delete pMySolver;
        delete pStimulus;
        VecDestroy(currentVoltage);
    }   
};
#endif //_TESTPARALLELMONODOMAINDG0ASSEMBLER_HPP_
