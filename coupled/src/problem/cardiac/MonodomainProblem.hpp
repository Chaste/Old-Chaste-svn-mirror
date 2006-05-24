#ifndef MONODOMAINPROBLEM_HPP_
#define MONODOMAINPROBLEM_HPP_

//#include <iostream>
// For mkdir()
#include <sys/stat.h>
#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "BoundaryConditionsContainer.hpp"
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ParallelColumnDataWriter.hpp"
#include "MonodomainPde.hpp"
#include "AbstractCardiacCellFactory.hpp"

/**
 * Class which specifies and solves a monodomain problem.
 */

template<int SPACE_DIM>
class MonodomainProblem
{
private:
    std::string mMeshFilename;
    double mStartTime;
    double mEndTime;
    std::string  mOutputDirectory, mOutputFilenamePrefix;

    MonodomainPde<SPACE_DIM> *mpMonodomainPde;
    bool mDebugOn;
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
     
    MonodomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
    : mMeshFilename(""),   // i.e. undefined
      mEndTime(1000),   // 1,000 ms = 1 second
      mOutputDirectory(""),   // i.e. undefined
      mOutputFilenamePrefix(""),   // i.e. undefined
      mpMonodomainPde(NULL),
      mDebugOn(false)
    {
        mpCellFactory = pCellFactory;
                
        mStartTime   = 0.0;  // ms
        mPdeTimeStep = 0.01; // ms
    }

    /**
     * Destructor
     */
     
    ~MonodomainProblem()
    { 
        if (mpMonodomainPde != NULL)
        {
            delete mpMonodomainPde;
        }
    }
    
    
    void Initialise()
    {
        assert( mMeshFilename!="" );

        TrianglesMeshReader mesh_reader(mMeshFilename);
        mMesh.ConstructFromMeshReader(mesh_reader);
        
        mpCellFactory->SetMesh( &mMesh );
        
        mpMonodomainPde = new MonodomainPde<SPACE_DIM>( mpCellFactory, mStartTime, mPdeTimeStep);
    }
     
    /**
     * Solve the problem
     */
    void Solve()
    {
        assert( mpMonodomainPde != NULL );

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
        MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM> monodomain_assembler(&linear_solver);
        monodomain_assembler.SetMatrixIsConstant();
        
        // initial condition;   
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, mMesh.GetNumNodes() );
        VecSetFromOptions(initial_condition);
  
        double* p_initial_condition;
        VecGetArray(initial_condition, &p_initial_condition); 
        
        VecGetOwnershipRange(initial_condition, &mLo, &mHi);
        
        // Set a constant initial voltage throughout the mMesh
        for (int global_index=mLo; global_index<mHi; global_index++)
        {
            int local_index = global_index - mLo;
            p_initial_condition[local_index] = mpMonodomainPde->GetCardiacCell(global_index)->GetVoltage();
        }
        VecRestoreArray(initial_condition, &p_initial_condition);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);
    
        /*
         *  Write data to a file <mOutputFilenamePrefix>_xx.dat, 'xx' refers to 
         *  'xx'th time step using ColumnDataWriter 
         */         
        
        ParallelColumnDataWriter *p_test_writer = NULL;
       
        int time_var_id = 0;
        int voltage_var_id = 0;
        if (mOutputFilenamePrefix.length() > 0)
        {        
            mkdir(mOutputDirectory.c_str(), 0777);
                 
            p_test_writer = new ParallelColumnDataWriter(mOutputDirectory,mOutputFilenamePrefix);

            p_test_writer->DefineFixedDimension("Node", "dimensionless", mMesh.GetNumNodes() );
            time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
        
            voltage_var_id = p_test_writer->DefineVariable("V","mV");
            p_test_writer->EndDefineMode();
        }
        else
        {
            throw ("mOutputFilenamePrefix should not be the empty string");
        }
        double current_time = mStartTime;        
        int big_steps = 0;
        
        p_test_writer->PutVariable(time_var_id, current_time); 
        p_test_writer->PutVector(voltage_var_id, initial_condition);
                        
        
        while( current_time < mEndTime )
        {
            monodomain_assembler.SetTimes(current_time, current_time+mPdeTimeStep, mPdeTimeStep);
            monodomain_assembler.SetInitialCondition( initial_condition );
            
            mVoltage = monodomain_assembler.Solve(mMesh, mpMonodomainPde, bcc);
            
            // Free old initial condition
            VecDestroy(initial_condition);

            // Initial condition for next loop is current solution
            initial_condition = mVoltage;
            
            // Writing data out to the file <mOutputFilenamePrefix>.dat                 
            p_test_writer->AdvanceAlongUnlimitedDimension(); //creates a new file
            p_test_writer->PutVariable(time_var_id, current_time); 
            p_test_writer->PutVector(voltage_var_id, mVoltage);


            /* //WARNING: won't run in parallel                   
            if (mDebugOn==true)
            {
                double* p_voltage;
                VecGetArray(mVoltage, &p_voltage);                    

                double max_voltage=p_voltage[0];
                int max_index=0;

                for(int j=0; j<mMesh.GetNumNodes(); j++) 
                {
                    if (p_voltage[j]>max_voltage)
                    {
                        max_voltage=p_voltage[j];
                        max_index=j;
                    }
                } 
                std::cout<<"At time "<< current_time <<" max voltage is "<<max_voltage<<" at "<<max_index<<"\n";          
                VecRestoreArray(mVoltage, &p_voltage); 
            }
            */
            
            mpMonodomainPde->ResetAsUnsolvedOdeSystem();

            current_time += mPdeTimeStep;
                
            big_steps++;
        }


        // close the file that stores voltage values            
        p_test_writer->Close();
        delete p_test_writer;
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
    
    MonodomainPde<SPACE_DIM> * GetMonodomainPde() 
    {
        return mpMonodomainPde;  
    }
};

#endif /*MONODOMAINPROBLEM_HPP_*/
