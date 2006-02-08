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
#include <cmath>
#include "EulerIvpOdeSolver.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestMonodomainFitzHughNagumoWithDg0Assembler : public CxxTest::TestSuite 
{   
public:
    
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
        EulerIvpOdeSolver mySolver;
        MonodomainPdeFitzHughNagumo<1> monodomain_pde(mesh.GetNumNodes(), &mySolver, tStart, tBigStep, tSmallStep);
        
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
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus);
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
        Vec initial_condition; 
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
      
        VecSetSizes(initial_condition, PETSC_DECIDE, mesh.GetNumNodes() );
        //VecSetType(initial_condition, VECSEQ);
        VecSetFromOptions(initial_condition);
  
        double* initial_condition_array;
        int lo, hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        int ierr = VecGetArray(initial_condition, &initial_condition_array); 
        
        // initial voltage condition of a constant everywhere on the mesh
        // ALREADY CHECKED
        for (int global_index=lo; global_index<hi; global_index++)
        {
            // change intiial conditions to exp(-x^2/10)
            double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
            //std::cout << i << " " << x << std::endl;
            initial_condition_array[global_index-lo] = exp(-(x*x)/10);
        }
        VecRestoreArray(initial_condition, &initial_condition_array);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);

        /*
         * Write data to a file NewMonodomainFHN_1d_xx.dat, 'xx' refers to nth time step
         *  using ColumnDataWriter 
         */                                                                            
        
        // Uncomment all column writer related lines to write data
        ColumnDataWriter *writer;
        writer = new ColumnDataWriter("testoutput","NewMonodomainFHN_1d");
       
        int time_var_id = 0;
        int voltage_var_id = 0;

        writer->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
        time_var_id = writer->DefineUnlimitedDimension("Time","msecs");

        voltage_var_id = writer->DefineVariable("V","mV");
        writer->EndDefineMode();
           
        double tCurrent = tStart;  
        int counter = 0;      
        Vec currentVoltage;
        double *currentVoltageArray;
        
        while( tCurrent < tFinal )
        {

            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( initial_condition );
            try {
                currentVoltage = monodomainAssembler.Solve(mesh, &monodomain_pde, bcc, &linearSolver);
             } catch (Exception e) {
                TS_TRACE(e.GetMessage());
                TS_ASSERT(0);
            }
            
            // Free old initial condition
            VecDestroy(initial_condition);
            // Initial condition for next loop is current solution
            initial_condition = currentVoltage;

            // Writing data out to the file NewMonodomainLR91_1d.dat
            /** \todo 
             * Here's a simple hack to make sure that the
             * data writers don't fall over in parallel.
             * It does *not* make them work!
             */
            if (counter % 20 == 0 && lo==0)    
            {
                int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
                writer->PutVariable(time_var_id, tCurrent); 
                // TS_TRACE("Put out voltage");
                for(int j=0; j</*mesh.GetNumNodes()*/ hi-lo; j++) 
                {
                    writer->PutVariable(voltage_var_id, currentVoltageArray[j], j);    
                    //std::cout << currentVoltageArray[j] << "\n" ;
                }
                
                VecRestoreArray(currentVoltage, &currentVoltageArray); 
                writer->AdvanceAlongUnlimitedDimension();
            } //end if currentTime
            
            monodomain_pde.ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
            counter++;
        }
        
        // close the file that stores voltage values
        writer->Close();
        delete writer;
        
        // test whether voltage is in correct range
        ierr = VecGetArray(currentVoltage, &currentVoltageArray);
        
        for(int local_index=0; local_index<hi-lo; local_index++)
        {
            TS_ASSERT_LESS_THAN_EQUALS(currentVoltageArray[local_index], 50);
            TS_ASSERT_LESS_THAN_EQUALS(-100, currentVoltageArray[local_index]);
        }
        VecRestoreArray(currentVoltage, &currentVoltageArray);

        VecDestroy(currentVoltage);
    }

}; 
#endif //_TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_HPP_
