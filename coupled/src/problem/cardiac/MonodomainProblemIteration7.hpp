#ifndef _MONODOMAINPROBLEMITERATION7_HPP_
#define _MONODOMAINPROBLEMITERATION7_HPP_

#include <iostream>

#include "MonodomainProblemIteration7.hpp"

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

#include "AbstractMonodomainProblemStimulus.hpp"

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
    bool mSequential; 
    bool mDebugOn;
    
    double mEndTime;
    AbstractMonodomainProblemStimulus<SPACE_DIM> *mpStimulus;
    std::string mMeshFilename, mOutputDirectory, mOutputFilenamePrefix;
    double mPdeTimeStep;  //aka big_timestep
    double mOdeTimeStep;  //aka small_timestep or ickle_timestep(jameso)
    // mOdeTimeStep will have to b removed when ticket 69 has been dealt with.
    
public:
    Vec mCurrentVoltage; // Current solution
    int mLo, mHi;
    MonodomainPde<SPACE_DIM> *mMonodomainPde;
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> mMesh;
    
    
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
     
//    MonodomainProblemIteration7(const std::string &rMeshFilename,
//                      const double &rEndTime,
//                      const std::string &rOutputDirectory,
//                      const std::string &rOutputFilenamePrefix,
//                      AbstractMonodomainProblemStimulus<SPACE_DIM> *rStimulus,
//                      const bool& rContainsInternalFaces = true,
//                      const bool& rDebug = false)
    MonodomainProblemIteration7()
    : mMeshFilename(""),   // i.e. undefined
      mEndTime(1000),   // 1,000 ms = 1 second
      mOutputDirectory(""),   // i.e. undefined
      mOutputFilenamePrefix(""),   // i.e. undefined
      mpStimulus(NULL),   // i.e. none
      mMonodomainPde(NULL),
      mDebugOn(false)
    {
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        mSequential = (num_procs == 1);
        mPdeTimeStep = 0.01; // ms
        mOdeTimeStep = 0.01; // ms
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
    
     
    /**
     * Solve the problem
     */
    void Solve(const double& rDiffusionCoefficient = 0.0005)
    {
        try
        {
            double start_time = 0.0;
        
            //double big_time_step = time_step; 
            //double small_time_step = time_step/2.0;
        
            // Read mMesh
            TrianglesMeshReader mesh_reader(mMeshFilename);
            mMesh.ConstructFromMeshReader(mesh_reader);
        
            // Instantiate PDE object
            MockEulerIvpOdeSolver ode_solver;
            mMonodomainPde = new MonodomainPde<SPACE_DIM>(mMesh.GetNumNodes(), &ode_solver, start_time, mPdeTimeStep, mOdeTimeStep);

            // Set the diffusion coefficient
        
            mMonodomainPde->SetDiffusionCoefficient(rDiffusionCoefficient);

            // Add initial stim       
            mpStimulus->Apply(mMonodomainPde, &mMesh);
        
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
                monodomain_assembler.SetTimes(current_time, current_time+mPdeTimeStep, mPdeTimeStep);
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
                    VecRestoreArray(mCurrentVoltage, &p_current_voltage); 
                }
                
                if (mDebugOn==true)
                {
                    odeVariablesType ode_vars;
                    int node_number;
                    
                    node_number=37876;
                    ode_vars = mMonodomainPde->GetOdeVarsAtNode(node_number);
                    std::cout<<"At time "<<current_time<<" node "<<node_number<<":\t"; 
                    for (unsigned i=0; i<ode_vars.size(); i++){
                        std::cout<<"("<<i<<") "<<ode_vars[i]<<"\t";   
                    }
                    std::cout<<std::endl;
                    
                    node_number=37877;
                    ode_vars = mMonodomainPde->GetOdeVarsAtNode(node_number);
                    std::cout<<"At time "<<current_time<<" node "<<node_number<<":\t"; 
                    for (unsigned i=0; i<ode_vars.size(); i++){
                        std::cout<<"("<<i<<") "<<ode_vars[i]<<"\t";   
                    }
                    std::cout<<std::endl;
                }
                
                mMonodomainPde->ResetAsUnsolvedOdeSystem();
                current_time += mPdeTimeStep;
                    
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
    
    void SetOdeTimeStep(double odeTimeStep)
    {
        assert(0.0 < odeTimeStep && odeTimeStep <= mPdeTimeStep);
        mOdeTimeStep=odeTimeStep;
    }

    void SetPdeTimeStep(double pdeTimeStep)
    {
        assert(0.0 < pdeTimeStep &&  mOdeTimeStep <= pdeTimeStep);
        mPdeTimeStep=pdeTimeStep;
    }

    void SetTimeSteps(double odeTimeStep, double pdeTimeStep)
    {
        assert(0.0 < odeTimeStep &&  0.0 < pdeTimeStep &&  odeTimeStep <= pdeTimeStep);
        mPdeTimeStep=pdeTimeStep;
        mOdeTimeStep=odeTimeStep;
    }
    
 
    double GetOdeTimeStep()
    {
        return mOdeTimeStep;   
    }
       
    double GetPdeTimeStep()
    {
        return mPdeTimeStep;   
    }

    void SetMeshFilename(const std::string &rMeshFilename)
    {
        mMeshFilename = rMeshFilename;
    }
    
    void SetEndTime(const double &rEndTime)
    {
        mEndTime = rEndTime;
    }
        
    void SetOutputDirectory(const std::string &rOutputDirectory)
    {
        mOutputDirectory = rOutputDirectory;
    }
        
    void SetOutputFilenamePrefix(const std::string &rOutputFilenamePrefix)
    {
        mOutputFilenamePrefix = rOutputFilenamePrefix;
    }
        
    void SetStimulus(AbstractMonodomainProblemStimulus<SPACE_DIM> *rStimulus)
    {
        mpStimulus = rStimulus;
    }
};


#endif //_MONODOMAINPROBLEMITERATION7_HPP_
