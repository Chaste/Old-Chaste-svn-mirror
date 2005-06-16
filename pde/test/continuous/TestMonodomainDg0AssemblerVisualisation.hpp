#ifndef _TESTMONODOMAINDG0ASSEMBLERVISUALISATION_HPP_
#define _TESTMONODOMAINDG0ASSEMBLERVISUALISATION_HPP_

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
#include "MonodomainPde.hpp"
#include "MonodomainDg0Assembler.hpp"
#include "ColumnDataWriter.hpp"
#include <cmath>
#include "MatlabVisualizer.cpp"

 
class TestMonodomainDg0AssemblerVisualisation : public CxxTest::TestSuite 
{   
public:
    void setUp()
    {
        int FakeArgc=0;
        char *FakeArgv0="testrunner";
        char **FakeArgv=&FakeArgv0;
        
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
    }   

    void TestMonodomainDg01DVisual()
    {  
        
        double tStart = 0; 
        double tFinal = 10;//400;//0.1;
        
        // use big time step (the pde timestep) is the same as the small time step (the ode timestep)
        double tBigStep = 0.01; 
        double tSmallStep = 0.01;
        
        // Create mesh from mesh reader 
		// TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
        
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        MonodomainPde<1> monodomain_pde(mesh.GetNumNodes(), pMySolver, tStart, tBigStep, tSmallStep);
        
        // sets Luo Rudy system with initial conditions passed on
        double voltage = -9999; // This voltage will be ignored
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;  
        double d = 0.003;
        double f = 1.0;
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
        Vec currentVoltage;
        VecCreate(PETSC_COMM_WORLD, &currentVoltage);
        VecSetSizes(currentVoltage, PETSC_DECIDE, mesh.GetNumNodes() );
        //VecSetType(initialCondition, VECSEQ);
        VecSetFromOptions(currentVoltage);
  
        double* currentVoltageArray;
        int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        // initial voltage condition of a constant everywhere on the mesh
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            currentVoltageArray[i] = -84.5;
        }
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(currentVoltage);
        VecAssemblyEnd(currentVoltage);

              
        /*
         * Write data to a file NewMonodomainLR91_1d_xx.dat, 'xx' refers to nth time step
         *  using ColumnDataWriter 
         */                                                                            
           
        
		//uncomment all column writer related lines to write data
        ColumnDataWriter *mpTestWriter;
        mpTestWriter = new ColumnDataWriter("data","NewMonodomainLR91_1d");
       
        int time_var_id = 0;
        int voltage_var_id = 0;

        mpTestWriter->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
        mpTestWriter->DefineUnlimitedDimension("Time","msecs");

        time_var_id = mpTestWriter->DefineVariable("Time","msecs");
        voltage_var_id = mpTestWriter->DefineVariable("V","mV");
        mpTestWriter->EndDefineMode();
           
        double tCurrent = tStart;  
        int counter;      
        while( tCurrent < tFinal )
        { 
            // std::cout << "t = " << tCurrent << "\n" << std::flush;

            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( currentVoltage );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            //std::cout << "Current time =  "<< tCurrent << "\n"; 
            // Writing data out to the file NewMonodomainLR91_1d.dat
     		
			if (counter % 20 == 0)    
            {
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
			} //end if currentTime
	        monodomain_pde.ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
			counter++;
        }
		
        // close the file that stores voltage values
        mpTestWriter->Close();
		delete mpTestWriter;

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
					//std::cout<< "gating variable is "<< odeVars[j] << sdt::endl;
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
                //test that Ca concentration is always great than 
                else if (j==3)
                {
                	TS_ASSERT_LESS_THAN_EQUALS ( -odeVars[3], 0.0);
                }
            }
        }
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
		VecDestroy(currentVoltage);
		
		// Visualise using Matlab
		AbstractVisualizer<1> *pViewer;
		TS_ASSERT_THROWS_NOTHING(
								 pViewer=new MatlabVisualizer<1>(
										"data/NewMonodomainLR91_1d"));
		try {
			pViewer->CreateFilesForVisualization();
		} catch (Exception e) {
			TS_TRACE(e.getMessage());
		}
		delete pViewer;
		
		// Tidy up
		delete pMySolver;
		delete pStimulus;
	}
    
    void NotestMonodomainDg02DVISUAL( void )
    {   
        double tStart = 0.0; 
        double tFinal = 10.0;
        
        // use big time step (the pde timestep) is the same as the small time step (the ode timestep)
        double tBigStep = 0.01; 
        double tSmallStep  = 0.01;
        
        // read mesh on [0,1]x[0,1]
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        MonodomainPde<2> monodomain_pde(mesh.GetNumNodes(), pMySolver, tStart, tBigStep, tSmallStep);
                
        
        // sets Luo Rudy system with initial conditions passed on
        double voltage = -9999; // This voltage will be ignored
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;
        double d = 0.003;
        double f = 1.0;
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
        MonodomainDg0Assembler<2,2> monodomainAssembler;
        
        // initial condition;   
        Vec currentVoltage;
        VecCreate(PETSC_COMM_WORLD, &currentVoltage);
        VecSetSizes(currentVoltage, PETSC_DECIDE, mesh.GetNumNodes() );
        //VecSetType(initialCondition, VECSEQ);
        VecSetFromOptions(currentVoltage);
  
        double* currentVoltageArray;
        int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        // initial voltage condition of a constant everywhere on the mesh
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            currentVoltageArray[i] = -84.5;
        }
     
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(currentVoltage);
        VecAssemblyEnd(currentVoltage);

              
        /*
        * Write data to a file NewMonodomainLR91_2d_xx.dat, 'xx' refers to nth time step
        *  using ColumnDataWriter 
        */                                                           
                 
     
        
//        // uncomment all column writer related lines to write data (and the line further below)         
        ColumnDataWriter *mpTestWriter;
        mpTestWriter = new ColumnDataWriter("data","NewMonodomainLR91_2d");
       
        int time_var_id = 0;
        int voltage_var_id = 0;

        mpTestWriter->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
        mpTestWriter->DefineUnlimitedDimension("Time","msecs");

        time_var_id = mpTestWriter->DefineVariable("Time","msecs");
        voltage_var_id = mpTestWriter->DefineVariable("V","mV");
        mpTestWriter->EndDefineMode();
//           
        double tCurrent = tStart;      
        int counter = 0;  
        while( tCurrent < tFinal )
        {
            // std::cout << "t = " << tCurrent << "\n" << std::flush;

            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( currentVoltage );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            
            // Writing data out to the file NewMonodomainLR91_2d.dat
            if (counter % 20 == 0.0)    
	        {
	            int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
	             
	              
	            mpTestWriter->PutVariable(time_var_id, tCurrent); 
	            for(int j=0; j<mesh.GetNumNodes(); j++) 
	            {
	                //uncomment the line below to write data (and the line further below)
	                mpTestWriter->PutVariable(voltage_var_id, currentVoltageArray[j], j);    
	            }
	  
	            VecRestoreArray(currentVoltage, &currentVoltageArray); 
	//             //uncomment the line below to write data
	            mpTestWriter->AdvanceAlongUnlimitedDimension();
            }
            counter ++;
            monodomain_pde.ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
        }

//        // close the file that stores voltage values
        mpTestWriter->Close();
    
            
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
 


    }   

    //uncommen these lines to visualize data using Matlab.     
 	void NotestVisualization2DNewMonodomainLR91Visual(void)
	{
		AbstractVisualizer<2> *pViewer;
		TS_ASSERT_THROWS_NOTHING(
		                  pViewer=new MatlabVisualizer<2>(
		                  "data/NewMonodomainLR91_2d"));
		TS_ASSERT_THROWS_NOTHING(
					pViewer->CreateFilesForVisualization());
		delete pViewer;

	}
	
    void NOtestMonodomainDg03DVisual( void )
    {   
        double tStart = 0; 
        double tFinal = 0.1;
         
        
        double tBigStep = 0.01; 
        double tSmallStep  = 0.01;
        
        // read 3d mesh
        // TrianglesMeshReader mesh_reader("pdes/tests/meshdata/cylinder_with_hole_840_elements");
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/slab_138_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        
        // Instantiate PDE object
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        MonodomainPde<3> monodomain_pde(mesh.GetNumNodes(), pMySolver, tStart, tBigStep, tSmallStep);
                
        
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
        MonodomainDg0Assembler<3,3> monodomainAssembler;
        
        // initial condition;   
        Vec currentVoltage;
        VecCreate(PETSC_COMM_WORLD, &currentVoltage);
        VecSetSizes(currentVoltage, PETSC_DECIDE, mesh.GetNumNodes() );
        //VecSetType(initialCondition, VECSEQ);
        VecSetFromOptions(currentVoltage);
  
        double* currentVoltageArray;
        int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
        // initial voltage condition of a constant everywhere on the mesh
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            currentVoltageArray[i] = -84.5;
        }
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(currentVoltage);
        VecAssemblyEnd(currentVoltage);

              
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
            monodomainAssembler.SetInitialCondition( currentVoltage );
            
            currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            
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
           
    }   
};
#endif //_TESTMONODOMAINDG0ASSEMBLERVISUALISATION_HPP_
