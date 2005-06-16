#ifndef _TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_HPP_

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

#include "MonodomainPdeFitzHughNagumo.hpp"

#include "MonodomainDg0Assembler.hpp"
#include "ColumnDataWriter.hpp"
#include "math.h"
#include "MatlabVisualizer.cpp"

 
class TestMonodomainFitzHughNagumoWithDG0Assembler : public CxxTest::TestSuite 
{   
public:
    void setUp()
    {
        int FakeArgc=0;
        char *FakeArgv0="testrunner";
        char **FakeArgv=&FakeArgv0;
        
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
    }
    
    /** \todo Not yet fully working...
     */
    void TestVisualFHN1D()
    {  
        
        double tStart = 0.0; 
        double tFinal = 10;//400;//0.1;
        
        // use big time step (the pde timestep) is the same as the small time step (the ode timestep)
        double tBigStep = 0.01; 
        double tSmallStep  = 0.01;
        
        // Create mesh from mesh reader
        //TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
        TrianglesMeshReader mesh_reader("mesh/test/data/heart_FHN_mesh");
        //TrianglesMeshReader mesh_reader("pdes/tests/meshdata/trivial_1d_mesh");
  
        //TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
       	MonodomainPdeFitzHughNagumo<1> monodomain_pde(mesh.GetNumNodes(), pMySolver, tStart, tBigStep, tSmallStep);
        
        // sets FHN system with initial conditions passed on
        double voltage = -9999; // This voltage will be ignored
        double w = 0.0;         // recovery variable
        
        double magnitudeOfStimulus = 1.0;  
        double durationOfStimulus  = 0.5 ;
                  
        // bad 
        std::vector<double> initialConditions;
        initialConditions.push_back(voltage);
        initialConditions.push_back(w);

        // set this as the initial condition of the gating vars at each node in the mesh        
        monodomain_pde.SetUniversalInitialConditions(initialConditions);
        
        // add initial stim to node 0 only
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus);
        // NO stimulus applied
        //        monodomain_pde.SetStimulusFunctionAtNode(5, pStimulus);
                
        
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
        // ALREADY CHECKED
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
        	// change intiial conditions to exp(-x^2/10)
        	double x = mesh.GetNodeAt(i)->GetPoint()[0];
        	//std::cout << i << " " << x << std::endl;
            currentVoltageArray[i] = exp(-(x*x)/10);
        }
        VecRestoreArray(currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(currentVoltage);
        VecAssemblyEnd(currentVoltage);

        /*
         * Write data to a file NewMonodomainFHN_1d_xx.dat, 'xx' refers to nth time step
         *  using ColumnDataWriter 
         */                                                                            
        
		// Uncomment all column writer related lines to write data
		ColumnDataWriter *mpTestWriter;
        mpTestWriter = new ColumnDataWriter("testoutput","NewMonodomainFHN_1d");
       
        int time_var_id = 0;
        int voltage_var_id = 0;

        mpTestWriter->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
        mpTestWriter->DefineUnlimitedDimension("Time","msecs");

        time_var_id = mpTestWriter->DefineVariable("Time","msecs");
        voltage_var_id = mpTestWriter->DefineVariable("V","mV");
        mpTestWriter->EndDefineMode();
           
        double tCurrent = tStart;  
        int counter = 0;      
        
        while( tCurrent < tFinal )
        {

            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( currentVoltage );
            try {
            	currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
            } catch (Exception e) {
            	TS_TRACE(e.getMessage());
            }

            // Writing data out to the file NewMonodomainLR91_1d.dat
			if (counter % 20 == 0)    
			{
				int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
				mpTestWriter->PutVariable(time_var_id, tCurrent); 
				// TS_TRACE("Put out voltage");
				for(int j=0; j<mesh.GetNumNodes(); j++) 
				{
					mpTestWriter->PutVariable(voltage_var_id, currentVoltageArray[j], j);    
					//std::cout << currentVoltageArray[j] << "\n" ;
				}
				
				VecRestoreArray(currentVoltage, &currentVoltageArray); 
				mpTestWriter->AdvanceAlongUnlimitedDimension();
			} //end if currentTime
			
			
			monodomain_pde.ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
			counter++;
        }
		
        // close the file that stores voltage values
        mpTestWriter->Close();
		delete mpTestWriter;
		
//        // test whether voltages and gating variables are in correct ranges
//        ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
//        
//        for(int i=0; i<mesh.GetNumNodes(); i++)
//        {
//            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
//            double Ena   =  54.4;
//            double Ek    = -77.0;
//            double Istim = -80.0;
//            
//            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[i] , Ena +  30);
//            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[i] + (Ek-30), 0);
//                
//            std::vector<double> odeVars = monodomain_pde.GetOdeVarsAtNode(i);           
//            for(int j=0; j<8; j++)
//            {
//                // if not voltage or calcium ion conc, test whether between 0 and 1 
//                if((j!=4) && (j!=3))
//                {
//                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
//                   //std::cout<< "gating variable is "<< odeVars[j] << sdt::endl;                   
//                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
//                }
//                //test that Ca concentration is always great than 
//                else if (j==3)
//                {
//                	TS_ASSERT_LESS_THAN_EQUALS ( -odeVars[3], 0.0);
//                }
//            }
//        }
//        VecRestoreArray(currentVoltage, &currentVoltageArray);      

	VecDestroy(currentVoltage);

	

	// Generate output files for visualising in Matlab
	AbstractVisualizer<1> *pViewer;
	TS_ASSERT_THROWS_NOTHING(
				 pViewer=new MatlabVisualizer<1>(
		                  "testoutput/NewMonodomainFHN_1d"));
	//TS_ASSERT_THROWS_NOTHING(pViewer->CreateFilesForVisualization());	
	try {
	    pViewer->CreateFilesForVisualization();
	} catch (Exception e) {
	    TS_TRACE(e.getMessage());
	}
	delete pViewer;

    }
    

}; // _TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_

#endif


