#ifndef _TESTMONODOMAINDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINDG0ASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
#include <iostream>
#include <cmath>
#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"  
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ColumnDataWriter.hpp"

#include "MonodomainPde.hpp"
#include "FischerPde.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestMonodomainDg0Assembler : public CxxTest::TestSuite 
{   
public:
	
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
	
	void TestMonodomainDg0AssemblerWithFischer1DAgainstSimpleDg0Assembler()
	{
		
		double tStart = 0;
		double tFinal = 1;
		double tBigStep = 0.01;
		// Create mesh from mesh reader 
		TrianglesMeshReader mesh_reader("mesh/test/data/heart_FHN_mesh");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
        
		// Instantiate PDE object
		FischerPde<1> pde;
        
		// Boundary conditions: zero neumann on entire boundary (2 elements)
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		ConstBoundaryCondition<1>* pNeumannBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
		ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
		bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);
		iter++;
		bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);

		// Linear solver
		SimpleLinearSolver linearSolver;

		// Assembler
		MonodomainDg0Assembler<1,1> monodomain_assembler;
		SimpleDg0ParabolicAssembler<1,1> simple_assembler;
        
		// initial condition;   
		Vec initial_condition_1, initial_condition_2;
		initial_condition_1 = CreateInitialConditionVec(mesh.GetNumNodes());
		VecDuplicate(initial_condition_1, &initial_condition_2);
  
		double* init_array;
	    int lo, hi;
		VecGetOwnershipRange(initial_condition_1,&lo,&hi);
		int ierr = VecGetArray(initial_condition_1, &init_array); 
		for (int i=lo; i<hi; i++)
		{
			double x=mesh.GetNodeAt(i)->GetPoint()[0];
			init_array[i-lo] = exp(-(x*x)/100);
		}
		VecRestoreArray(initial_condition_1, &init_array);      
		VecAssemblyBegin(initial_condition_1);
		VecAssemblyEnd(initial_condition_1);
		//VecView(initial_condition_1, PETSC_VIEWER_STDOUT_WORLD);
		VecCopy(initial_condition_1, initial_condition_2); // Both assemblers use same initial cond'n

		// Vars to hold current solutions at each iteration
		Vec current_solution_1, current_solution_2;

		double tCurrent = tStart;
		while( tCurrent < tFinal )
		{
			monodomain_assembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
			simple_assembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
			
			monodomain_assembler.SetInitialCondition( initial_condition_1 );
			simple_assembler.SetInitialCondition( initial_condition_2 );

			//std::cout<<"monodomain solve\n";std::cerr.flush();std::cout.flush();
			current_solution_1 = monodomain_assembler.Solve(mesh, &pde, bcc, &linearSolver);
			//std::cout<<"simple solve\n";
			current_solution_2 = simple_assembler.Solve(mesh, &pde, bcc, &linearSolver);
			//std::cout<<" solved\n";
     
     		// Next iteration uses current solution as initial condition
     		VecDestroy(initial_condition_1); // Free old initial condition
     		VecDestroy(initial_condition_2);
     		initial_condition_1 = current_solution_1;
     		initial_condition_2 = current_solution_2;
     
			tCurrent += tBigStep;
		}
		
		// Compare the results
		double *res1, *res2;
		VecGetArray(current_solution_1, &res1);
		VecGetArray(current_solution_2, &res2);
		for (int i=lo; i<hi; i++)
		{
			TS_ASSERT_DELTA(res1[i-lo], res2[i-lo], 1e-3);
		}
		VecRestoreArray(current_solution_1, &res1);
		VecRestoreArray(current_solution_2, &res2);
		VecDestroy(current_solution_1);
		VecDestroy(current_solution_2);
	}
    
	void TestMonodomainDg01D()
	{
		double tStart = 0; 
		double tFinal = 5;
        
		// use big time step (the pde timestep) is the same as the small time step (the ode timestep)
		double tBigStep = 0.01; 
		double tSmallStep  = 0.01;
        
		// Create mesh from mesh reader 
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
        
		// Instantiate PDE object
		EulerIvpOdeSolver mySolver;
		MonodomainPde<1> monodomain_pde(mesh.GetNumNodes(), &mySolver, tStart, tBigStep, tSmallStep);
        
		// sets Luo Rudy system with initial conditions passed on
		double voltage = -9999; // This voltage will be ignored
		double m = 0.0017;
		double h = 0.9833;
		double j = 0.9895;  
		double d = 0.003;
		double f = 1;
		double x = 0.0056;
		double caI = 0.0002;
		double magnitudeOfStimulus = -600.0;  
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
		InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus);
		monodomain_pde.SetStimulusFunctionAtNode(0, &stimulus);
                
        
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
		MonodomainDg0Assembler<1,1> monodomainAssembler;
        
		// initial condition;   
		Vec initial_condition;
		VecCreate(PETSC_COMM_WORLD, &initial_condition);
		VecSetSizes(initial_condition, PETSC_DECIDE, mesh.GetNumNodes() );
		//VecSetType(initial_condition, VECSEQ);
		VecSetFromOptions(initial_condition);
  
		double* initial_condition_array;
		int ierr = VecGetArray(initial_condition, &initial_condition_array); 
        
		// initial voltage condition of a constant everywhere on the mesh
		for(int i=0; i<mesh.GetNumNodes(); i++)
		{
			initial_condition_array[i] = -84.5;
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		VecAssemblyBegin(initial_condition);
		VecAssemblyEnd(initial_condition);

              
		/*
		 * Write data to a file NewMonodomainLR91_1d_xx.dat, 'xx' refers to nth time step
		 *  using ColumnDataWriter 
		 */                                                                            
           
        
		// uncomment all column writer related lines to write data (and the line further below)         
		ColumnDataWriter *mpTestWriter;
		mpTestWriter = new ColumnDataWriter("testoutput","NewMonodomainLR91_1d");

		int time_var_id = 0;
		int voltage_var_id = 0;

		mpTestWriter->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
		mpTestWriter->DefineUnlimitedDimension("Time","msecs");

		time_var_id = mpTestWriter->DefineVariable("Time","msecs");
		voltage_var_id = mpTestWriter->DefineVariable("V","mV");
		mpTestWriter->EndDefineMode();
        
        Vec currentVoltage;  // Current solution
        double *currentVoltageArray;
        
		double tCurrent = tStart;        
		while( tCurrent < tFinal )
		{ 
            // std::cout << "t = " << tCurrent << "\n" << std::flush;

            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( initial_condition );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            
            // Free old initial condition
            VecDestroy(initial_condition);
            // Initial condition for next loop is current solution
            initial_condition = currentVoltage;
            
            // Writing data out to the file NewMonodomainLR91_1d.dat
         
            int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
              
            mpTestWriter->PutVariable(time_var_id, tCurrent); 
            for(int j=0; j<mesh.GetNumNodes(); j++) 
            {
                //uncomment the line below to write data (and the line further below)
                mpTestWriter->PutVariable(voltage_var_id, currentVoltageArray[j], j);    
            }
  
            VecRestoreArray(currentVoltage, &currentVoltageArray); 
            //uncomment the line below to write data
            mpTestWriter->AdvanceAlongUnlimitedDimension();
     
            monodomain_pde.ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
        }

        // close the file that stores voltage values
        mpTestWriter->Close();

        // test whether voltages and gating variables are in correct ranges        
        ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[i] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[i] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_pde.GetOdeVarsAtNode(i);           
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
        EulerIvpOdeSolver mySolver;
        MonodomainPde<2> monodomain_pde(mesh.GetNumNodes(), &mySolver, tStart, tBigStep, tSmallStep);
                
        
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
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus);
        monodomain_pde.SetStimulusFunctionAtNode(0, &stimulus);
                         
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
        MonodomainDg0Assembler<2,2> monodomainAssembler;
        
        // initial condition;   
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, mesh.GetNumNodes() );
        //VecSetType(initial_condition, VECSEQ);
        VecSetFromOptions(initial_condition);
  
        double* initial_condition_array;
        int ierr = VecGetArray(initial_condition, &initial_condition_array); 
        
        // initial voltage condition of a constant everywhere on the mesh
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            initial_condition_array[i] = -84.5;
        }
     
        VecRestoreArray(initial_condition, &initial_condition_array);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);

              
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
        
        Vec currentVoltage; // Current solution
        
        double tCurrent = tStart;        
        while( tCurrent < tFinal )
        {
            // std::cout << "t = " << tCurrent << "\n" << std::flush;

            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( initial_condition );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            
            // Free old initial condition
            VecDestroy(initial_condition);
            // Initial condition for next loop is current solution
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
        double *currentVoltageArray;
        ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[i] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[i] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_pde.GetOdeVarsAtNode(i);           
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
 		VecDestroy(currentVoltage);


    }   


    void testMonodomainDg03D( void )
    {   
        double tStart = 0; 
        double tFinal = 0.1;
         
        
        double tBigStep = 0.01; 
        double tSmallStep  = 0.01;
        
        // read 3d mesh
        // TrianglesMeshReader mesh_reader("mesh/test/data/cylinder_with_hole_840_elements");
        TrianglesMeshReader mesh_reader("mesh/test/data/slab_138_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        
        // Instantiate PDE object
        EulerIvpOdeSolver mySolver;
        MonodomainPde<3> monodomain_pde(mesh.GetNumNodes(), &mySolver, tStart, tBigStep, tSmallStep);
                
        
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
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus);
        monodomain_pde.SetStimulusFunctionAtNode(0, &stimulus);
 
                 
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
        MonodomainDg0Assembler<3,3> monodomainAssembler;
        
        // initial condition;   
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, mesh.GetNumNodes() );
        //VecSetType(initial_condition, VECSEQ);
        VecSetFromOptions(initial_condition);
  
        double* initial_condition_array;
        int ierr = VecGetArray(initial_condition, &initial_condition_array); 
        
        // initial voltage condition of a constant everywhere on the mesh
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            initial_condition_array[i] = -84.5;
        }
        VecRestoreArray(initial_condition, &initial_condition_array);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);

              
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
        
        Vec currentVoltage; // Current solution
        
        double tCurrent = tStart;        
        while( tCurrent < tFinal )
        {
            // std::cout << "t = " << tCurrent << "\n" << std::flush;

            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( initial_condition );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            
            // Free old initial condition
            VecDestroy(initial_condition);
            // Initial condition for next loop is current solution
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
        double *currentVoltageArray;
        ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[i] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[i] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_pde.GetOdeVarsAtNode(i);           
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
        VecDestroy(currentVoltage);
           
    }   
};
#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
