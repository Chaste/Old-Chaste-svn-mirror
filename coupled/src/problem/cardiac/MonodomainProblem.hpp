#ifndef _MONODOMAINPROBLEM_HPP_
#define _MONODOMAINPROBLEM_HPP_

#include <iostream>

//#include "petscvec.h"
//#include <vector>
//#include <cmath>
//#include <sys/stat.h>
//#include <sys/types.h>

#include "MonodomainProblem.hpp"

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
#include "MockEulerIvpOdeSolver.hpp"
//#include "FischerPde.hpp"

//#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "GenericCoupledProblemStimulus.hpp"

template<int SPACE_DIM>
class MonodomainProblem
{
private:
    double tFinal;
    GenericCoupledProblemStimulus<SPACE_DIM> *stimulus;

    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> mesh;

    std::string mMeshFilename, mOutputDirectory, mOutputFilenamePrefix;

public:
    Vec currentVoltage; // Current solution
    int lo, hi;
    MonodomainPde<SPACE_DIM> *monodomain_pde;
    
    /**
     * Constructor
     */
     
    MonodomainProblem(const std::string &rMeshFilename,
                      const double &rFinal,
                      const std::string &rOutputDirectory,
                      const std::string &rOutputFilenamePrefix,
                      GenericCoupledProblemStimulus<SPACE_DIM> *rStimulus)
    : mMeshFilename(rMeshFilename),
      tFinal(rFinal),
      mOutputDirectory(rOutputDirectory),
      mOutputFilenamePrefix(rOutputFilenamePrefix),
      stimulus(rStimulus),
      monodomain_pde(NULL)
    {
    }

    /**
     * Destructor
     */
     
    ~MonodomainProblem()
    { 
        if (monodomain_pde != NULL) delete monodomain_pde;
    }

    /**
     * Solve the problem
     */
    void Solve(void)
    {
        double tStart = 0.0;
    
        double tBigStep = 0.01; 
        double tSmallStep = 0.01;
    
        // Read mesh
        TrianglesMeshReader mesh_reader(mMeshFilename);
        mesh.ConstructFromMeshReader(mesh_reader);
    
        // Instantiate PDE object
        MockEulerIvpOdeSolver mySolver;
        monodomain_pde = new MonodomainPde<SPACE_DIM>(mesh.GetNumNodes(), &mySolver, tStart, tBigStep, tSmallStep);
    
        // Add initial stim       
        stimulus->Apply(monodomain_pde);    
    
        // Boundary conditions, zero neumann everywhere
        BoundaryConditionsContainer<SPACE_DIM,SPACE_DIM> bcc(1, mesh.GetNumNodes());
       
        // The 'typename' keyword is required otherwise the compiler complains
        // Not totally sure why!
        typename ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
        ConstBoundaryCondition<SPACE_DIM>* pNeumannBoundaryCondition = new ConstBoundaryCondition<SPACE_DIM>(0.0);
        
        while(iter < mesh.GetLastBoundaryElement())
        {
            bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);
            iter++;
        }
        
        // Linear solver
        SimpleLinearSolver linearSolver;
    
        // Assembler
        MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM> monodomainAssembler;
        
        // initial condition;   
            Vec initial_condition;
            VecCreate(PETSC_COMM_WORLD, &initial_condition);
            VecSetSizes(initial_condition, PETSC_DECIDE, mesh.GetNumNodes() );
            VecSetFromOptions(initial_condition);
      
            double* initial_condition_array;
            VecGetArray(initial_condition, &initial_condition_array); 
            
            VecGetOwnershipRange(initial_condition, &lo, &hi);
            
            // initial voltage condition of a constant everywhere on the mesh
        for(int global_index=lo; global_index<hi; global_index++)
        {
            initial_condition_array[global_index-lo] = -84.5;
        }
        VecRestoreArray(initial_condition, &initial_condition_array);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);
    
        /*
         *  Write data to a file NewMonodomainLR91_3d_xx.dat, 'xx' refers to nth time step
         *  using ColumnDataWriter 
         */         
        
        mkdir(mOutputDirectory.c_str(), 0777);
                 
        ColumnDataWriter *mpTestWriter;
        mpTestWriter = new ColumnDataWriter(mOutputDirectory,mOutputFilenamePrefix);
       
        int time_var_id = 0;
        int voltage_var_id = 0;
    
        mpTestWriter->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
        time_var_id = mpTestWriter->DefineUnlimitedDimension("Time","msecs");
    
        voltage_var_id = mpTestWriter->DefineVariable("V","mV");
        mpTestWriter->EndDefineMode();
        
        double* currentVoltageArray;
        
        double tCurrent = tStart;        

        int big_steps = 0;
        
        while( tCurrent < tFinal )
        {
            monodomainAssembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            monodomainAssembler.SetInitialCondition( initial_condition );
            
            currentVoltage = monodomainAssembler.Solve(mesh, monodomain_pde, bcc, &linearSolver);
            
            // Free old initial condition
            VecDestroy(initial_condition);
            // Initial condition for next loop is current solution
            initial_condition = currentVoltage;
            
            // Writing data out to the file <mOutputFilenamePrefix>.dat
             
            VecGetArray(currentVoltage, &currentVoltageArray);

            mpTestWriter->PutVariable(time_var_id, tCurrent); 
            for(int j=0; j<mesh.GetNumNodes(); j++) 
            {
                mpTestWriter->PutVariable(voltage_var_id, currentVoltageArray[j], j);    
            }
  
            VecRestoreArray(currentVoltage, &currentVoltageArray); 
             mpTestWriter->AdvanceAlongUnlimitedDimension();
     
            monodomain_pde->ResetAsUnsolvedOdeSystem();
            tCurrent += tBigStep;
                
            big_steps++;
        }

        TS_ASSERT_EQUALS(mySolver.GetCallCount(), (hi-lo)*big_steps);

            // close the file that stores voltage values
        mpTestWriter->Close();
        delete mpTestWriter;
    }
};
#endif //_MONODOMAINPROBLEM_HPP_
