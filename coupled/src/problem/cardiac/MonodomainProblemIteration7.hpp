#ifndef MONODOMAINPROBLEMITERATION7_HPP_
#define MONODOMAINPROBLEMITERATION7_HPP_

//#include <iostream>

// For mkdir()
#include <sys/stat.h>
#include <sys/types.h>

#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ColumnDataWriter.hpp"

#include "MonodomainPdeIteration7.hpp"
#include "MockEulerIvpOdeSolver.hpp"

#include "AbstractCardiacCellFactory.hpp"

/**
 * Class which specifies and solves a monodomain problem.
 */

template<int SPACE_DIM>
class MonodomainProblemIteration7
{
private:
    /**
     * Flag that is true when running on one processor
     */
    std::string mMeshFilename;
    double mStartTime;
    double mEndTime;
    std::string  mOutputDirectory, mOutputFilenamePrefix;

    MonodomainPdeIteration7<SPACE_DIM> *mMonodomainPde;
    bool mDebugOn;
    bool mSequential; 
    double mPdeTimeStep;  //aka big_timestep

    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    
    Vec mVoltage; // Current solution
    int mLo, mHi;

   
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> mMesh;
    
public:
    
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
     
    MonodomainProblemIteration7(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
    : mMeshFilename(""),   // i.e. undefined
      mEndTime(1000),   // 1,000 ms = 1 second
      mOutputDirectory(""),   // i.e. undefined
      mOutputFilenamePrefix(""),   // i.e. undefined
      mMonodomainPde(NULL),
      mDebugOn(false)
    {
        mpCellFactory = pCellFactory;
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        mSequential = (num_procs == 1);
        
        mStartTime   = 0.0;  // ms
        mPdeTimeStep = 0.01; // ms
    }

    /**
     * Destructor
     */
     
    ~MonodomainProblemIteration7()
    { 
        if (mMonodomainPde != NULL)
        {
            delete mMonodomainPde;
        }
    }
    
    
    void Initialise()
    {
        assert( mMeshFilename!="" );

        TrianglesMeshReader mesh_reader(mMeshFilename);
        mMesh.ConstructFromMeshReader(mesh_reader);
        
        mpCellFactory->SetMesh( &mMesh );
        
        mMonodomainPde = new MonodomainPdeIteration7<SPACE_DIM>( mpCellFactory, mStartTime, mPdeTimeStep);
    }
     
    /**
     * Solve the problem
     */
    void Solve(const double& rDiffusionCoefficient = 0.0005)
    {
        assert( mMonodomainPde != NULL );
                
        try
        {
            // Set the diffusion coefficient
            mMonodomainPde->SetDiffusionCoefficient(rDiffusionCoefficient);

            // Boundary conditions, zero neumann everywhere
            BoundaryConditionsContainer<SPACE_DIM,SPACE_DIM> bcc(1, mMesh.GetNumNodes());
           
            // The 'typename' keyword is required otherwise the compiler complains
            // Not totally sure why!
            typename ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>::BoundaryElementIterator iter = mMesh.GetBoundaryElementIteratorBegin();
            ConstBoundaryCondition<SPACE_DIM>* p_neumann_boundary_condition = new ConstBoundaryCondition<SPACE_DIM>(0.0);
            
            while(iter != mMesh.GetBoundaryElementIteratorEnd())
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
            
            double current_time = mStartTime;        
    
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
                monodomain_assembler.SetTimes(current_time, current_time+mPdeTimeStep, mPdeTimeStep);
                monodomain_assembler.SetInitialCondition( initial_condition );
                
                mVoltage = monodomain_assembler.Solve(mMesh, mMonodomainPde, bcc, &linear_solver);
                
                // Free old initial condition
                VecDestroy(initial_condition);

                // Initial condition for next loop is current solution
                initial_condition = mVoltage;
                
                // Writing data out to the file <mOutputFilenamePrefix>.dat                 
                if (mSequential)
                {
                    p_test_writer->AdvanceAlongUnlimitedDimension(); //creates a new file
                    
                    VecGetArray(mVoltage, &p_current_voltage);
        
                    p_test_writer->PutVariable(time_var_id, current_time); 
                    
                    for(int j=0; j<mMesh.GetNumNodes(); j++) 
                    {
                        p_test_writer->PutVariable(voltage_var_id, p_current_voltage[j], j);    
                    }
                   
                   if (mDebugOn==true){
                       double max_voltage=p_current_voltage[0];
                       int max_index=0;
                       for(int j=0; j<mMesh.GetNumNodes(); j++) 
                        {
                            if (p_current_voltage[j]>max_voltage){
                                max_voltage=p_current_voltage[j];
                                max_index=j;
                            }
                        } 
                        std::cout<<"At time "<< current_time <<" max voltage is "<<max_voltage<<" at "<<max_index<<"\n";          
                   }
                    VecRestoreArray(mVoltage, &p_current_voltage); 
                }
                
                mMonodomainPde->ResetAsUnsolvedOdeSystem();
                current_time += mPdeTimeStep;
                    
                big_steps++;
            }
    
      //      TS_ASSERT_EQUALS(ode_solver.GetCallCount(), (mHi-mLo)*big_steps);
    
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
    


    void SetStartTime(const double &rStartTime)
    {
        mStartTime = rStartTime;
    }

    void SetEndTime(const double &rEndTime)
    {
        mEndTime = rEndTime;
    }

    void SetPdeTimeStep(double pdeTimeStep)
    {
        assert(0.0 < pdeTimeStep);
        mPdeTimeStep=pdeTimeStep;
    }
    
    double GetPdeTimeStep()
    {
        return mPdeTimeStep;   
    }

    void SetMeshFilename(const std::string &rMeshFilename)
    {
        mMeshFilename = rMeshFilename;
    }
        
    void SetOutputDirectory(const std::string &rOutputDirectory)
    {
        mOutputDirectory = rOutputDirectory;
    }
        
    void SetOutputFilenamePrefix(const std::string &rOutputFilenamePrefix)
    {
        mOutputFilenamePrefix = rOutputFilenamePrefix;
    }
    
    AbstractCardiacCell* GetCardiacCell( int globalIndex )
    {
        return mMonodomainPde->GetCardiacCell(globalIndex);
    }
    
    void GetVoltageArray(double **pVoltageArray, int &lo, int &hi)
    {
        VecGetArray(mVoltage, pVoltageArray);
        lo=mLo;
        hi=mHi; 
    }
    
    void RestoreVoltageArray(double **pVoltageArray)
    {
       VecRestoreArray(mVoltage, pVoltageArray);      
       VecAssemblyBegin(mVoltage);
       VecAssemblyEnd(mVoltage);
       VecDestroy(mVoltage);
    }
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> & rGetMesh(){
        return mMesh;   
    }
};



#endif /*MONODOMAINPROBLEMITERATION7_HPP_*/
