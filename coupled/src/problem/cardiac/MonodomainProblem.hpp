#ifndef _MONODOMAINPROBLEM_HPP_
#define _MONODOMAINPROBLEM_HPP_

#include <iostream>

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

#include "MonodomainProblem.hpp"
#include "AbstractMonodomainProblemStimulus.hpp"

/**
 * Class which specifies and solves a monodomain problem.
 */

template<int SPACE_DIM>
class MonodomainProblem
{
private:
    /**
     * Flag that is true when running on one processor
     */
    bool mSequential; 
    
    double mEndTime;
    AbstractMonodomainProblemStimulus<SPACE_DIM> *mpStimulus;
    std::string mMeshFilename, mOutputDirectory, mOutputFilenamePrefix;
    bool mContainsInternalFaces;

public:
    Vec mCurrentVoltage; // Current solution
    int mLo, mHi;
    MonodomainPde<SPACE_DIM> *mMonodomainPde;
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> mMesh;
    
    double time_step;
    
    /**
     * Constructor
     * @param rMeshFilename Name of mesh used in simulation.  Note that the space
     *     step is measured in cm.
     * @param rEndTime Duration of simulation, in milliseconds.
     * @param rOutputDirectory Directory where voltage for each time step is written.
     * @param rOutputFilePrefix Filename prefix for above. "_nnnnnn.dat" is appended where nnnnnn is the time step.
     * @param rStimulus Object specifying the stimulus information.
     * @param rContainsInternalFaces Optional parameter specifying whether the mesh contains internal faces. Default is true.
     */
     
    MonodomainProblem(const std::string &rMeshFilename,
                      const double &rEndTime,
                      const std::string &rOutputDirectory,
                      const std::string &rOutputFilenamePrefix,
                      AbstractMonodomainProblemStimulus<SPACE_DIM> *rStimulus,
                      const bool& rContainsInternalFaces = true)
    : mMeshFilename(rMeshFilename),
      mEndTime(rEndTime),
      mOutputDirectory(rOutputDirectory),
      mOutputFilenamePrefix(rOutputFilenamePrefix),
      mpStimulus(rStimulus),
      mMonodomainPde(NULL),
      mContainsInternalFaces(rContainsInternalFaces)
    {
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        mSequential = (num_procs == 1);
        time_step = 0.01; // ms
    }

    /**
     * Destructor
     */
     
    ~MonodomainProblem()
    { 
        if (mMonodomainPde != NULL)
        {
            delete mMonodomainPde;
        }
    }
    
     
    /**
     * Solve the problem
     */
    void Solve(void)
    {
        try
        {
            double start_time = 0.0;
        
            double big_time_step = time_step; 
            double small_time_step = time_step;
        
            // Read mMesh
            TrianglesMeshReader mesh_reader(mMeshFilename, mContainsInternalFaces);
            mMesh.ConstructFromMeshReader(mesh_reader);
        
            // Instantiate PDE object
            MockEulerIvpOdeSolver ode_solver;
            mMonodomainPde = new MonodomainPde<SPACE_DIM>(mMesh.GetNumNodes(), &ode_solver, start_time, big_time_step, small_time_step);
        
            // Add initial stim       
            mpStimulus->Apply(mMonodomainPde, &mMesh);
        
            // Boundary conditions, zero neumann everywhere
            BoundaryConditionsContainer<SPACE_DIM,SPACE_DIM> bcc(1, mMesh.GetNumNodes());
           
            // The 'typename' keyword is required otherwise the compiler complains
            // Not totally sure why!
            typename ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>::BoundaryElementIterator iter = mMesh.GetFirstBoundaryElement();
            ConstBoundaryCondition<SPACE_DIM>* p_neumann_boundary_condition = new ConstBoundaryCondition<SPACE_DIM>(0.0);
            
            while(iter < mMesh.GetLastBoundaryElement())
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);
                iter++;
            }
            
            // Linear solver
            SimpleLinearSolver linear_solver;
        
            // Assembler
            MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM> monodomain_assembler;
            monodomain_assembler.SetMatrixIsConstant(&linear_solver);
            
            
            // initial condition;   
            Vec initial_condition;
            VecCreate(PETSC_COMM_WORLD, &initial_condition);
            VecSetSizes(initial_condition, PETSC_DECIDE, mMesh.GetNumNodes() );
            VecSetFromOptions(initial_condition);
      
            double* p_initial_condition;
            VecGetArray(initial_condition, &p_initial_condition); 
            
            VecGetOwnershipRange(initial_condition, &mLo, &mHi);
            
            // Set a constant initial voltage throughout the mMesh
            for(int local_index=0; local_index<mHi - mLo; local_index++)
            {            	
                p_initial_condition[local_index] = -84.5;
            }
            VecRestoreArray(initial_condition, &p_initial_condition);      
            VecAssemblyBegin(initial_condition);
            VecAssemblyEnd(initial_condition);
        
            /*
             *  Write data to a file <mOutputFilenamePrefix>_xx.dat, 'xx' refers to 
             *  'xx'th time step using ColumnDataWriter 
             */         
            
            ColumnDataWriter *p_test_writer;
           
            int time_var_id = 0;
            int voltage_var_id = 0;

            if (mSequential && mOutputFilenamePrefix.length() > 0)
            {        
                mkdir(mOutputDirectory.c_str(), 0777);
                     
                p_test_writer = new ColumnDataWriter(mOutputDirectory,mOutputFilenamePrefix);
    
                p_test_writer->DefineFixedDimension("Node", "dimensionless", mMesh.GetNumNodes() );
                time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
            
                voltage_var_id = p_test_writer->DefineVariable("V","mV");
                p_test_writer->EndDefineMode();
            }
            
            double* p_current_voltage;
            
            double current_time = start_time;        
    
            int big_steps = 0;
            
            
           	if (mSequential)
            {
                p_test_writer->PutVariable(time_var_id, current_time); 
                for(int j=0; j<mMesh.GetNumNodes(); j++) 
                {
                    p_test_writer->PutVariable(voltage_var_id, p_initial_condition[j], j);    
                }
            }            
            
            
            while( current_time < mEndTime )
            {
                monodomain_assembler.SetTimes(current_time, current_time+big_time_step, big_time_step);
                monodomain_assembler.SetInitialCondition( initial_condition );
                
                mCurrentVoltage = monodomain_assembler.Solve(mMesh, mMonodomainPde, bcc, &linear_solver);
                
                // Free old initial condition
                VecDestroy(initial_condition);

                // Initial condition for next loop is current solution
                initial_condition = mCurrentVoltage;
                
                // Writing data out to the file <mOutputFilenamePrefix>.dat                 
                if (mSequential)
                {
                    p_test_writer->AdvanceAlongUnlimitedDimension(); //creates a new file
                    
                    VecGetArray(mCurrentVoltage, &p_current_voltage);
        
                    p_test_writer->PutVariable(time_var_id, current_time); 
                    for(int j=0; j<mMesh.GetNumNodes(); j++) 
                    {
                        p_test_writer->PutVariable(voltage_var_id, p_current_voltage[j], j);    
                    }          
                    VecRestoreArray(mCurrentVoltage, &p_current_voltage); 
                }
         
                mMonodomainPde->ResetAsUnsolvedOdeSystem();
                current_time += big_time_step;
                    
                big_steps++;
            }
    
            TS_ASSERT_EQUALS(ode_solver.GetCallCount(), (mHi-mLo)*big_steps);
    
            // close the file that stores voltage values            
            if (mSequential && mOutputFilenamePrefix.length() > 0)
            {
                p_test_writer->Close();
            
                delete p_test_writer;
            }
        }
        catch (Exception &e)
        {
            std::cout<<e.GetMessage()<<std::endl;   
        }
    }
};
#endif //_MONODOMAINPROBLEM_HPP_
