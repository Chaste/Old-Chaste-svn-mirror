#ifndef BIDOMAINPROBLEM_HPP_
#define BIDOMAINPROBLEM_HPP_


// For mkdir()
#include <sys/stat.h>
#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "BoundaryConditionsContainer.hpp"
#include "BidomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ParallelColumnDataWriter.hpp"
#include "BidomainPde.hpp"
#include "AbstractCardiacCellFactory.hpp"


/**
 * Class which specifies and solves a bidomain problem.
 */
template<int SPACE_DIM>
class BidomainProblem
{
private:
    std::string mMeshFilename;
    double mStartTime;
    double mEndTime;
    std::string  mOutputDirectory, mOutputFilenamePrefix;

    BidomainPde<SPACE_DIM>* mpBidomainPde;
    bool mDebugOn;
    double mPdeTimeStep;

    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    
    Vec mVoltage; // Current solution
    int mLo, mHi; // Ownership range for the double size potential vector mVoltage

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
     
    BidomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
    : mMeshFilename(""),     // i.e. undefined
      mEndTime(1000),        // 1,000 ms = 1 second
      mOutputDirectory(""),  // i.e. undefined
      mOutputFilenamePrefix(""),   // i.e. undefined
      mpBidomainPde(NULL),
      mDebugOn(false)
    {
        mpCellFactory = pCellFactory;
                
        mStartTime   = 0.0;  // ms
        mPdeTimeStep = 0.01; // ms
        
        //initialise these to -1
        mLo = -1;
        mHi = -1; 
    }

    /**
     * Destructor
     */
    ~BidomainProblem()
    { 
        if (mpBidomainPde != NULL)
        {
            delete mpBidomainPde;
        }
    }
    
    
    void Initialise()
    {
        assert( mMeshFilename!="" );

        TrianglesMeshReader mesh_reader(mMeshFilename);
        mMesh.ConstructFromMeshReader(mesh_reader);
        
        mpCellFactory->SetMesh( &mMesh );
        
        mpBidomainPde = new BidomainPde<SPACE_DIM>( mpCellFactory, mStartTime, mPdeTimeStep );
    }
     

    /**
     * Solve the problem
     */
    void Solve()
    {
        assert( mpBidomainPde != NULL );
        assert( mStartTime < mEndTime );
        
        // Linear solver
        SimpleLinearSolver linear_solver;
    
        // Assembler
        BidomainDg0Assembler<SPACE_DIM,SPACE_DIM> bidomain_assembler(mpBidomainPde, &mMesh, &linear_solver);
        
        // initial condition;   
        Vec initial_condition;
        int lo, hi;
        mpBidomainPde->GetOwnershipRange(lo, hi);
        VecCreateMPI(PETSC_COMM_WORLD, 2*(hi-lo) , 2*mMesh.GetNumNodes(), &initial_condition);  
        double* p_initial_condition;
        VecGetArray(initial_condition, &p_initial_condition); 

        VecGetOwnershipRange(initial_condition, &mLo, &mHi);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            p_initial_condition[2*local_index  ] = mpBidomainPde->GetCardiacCell(global_index)->GetVoltage();
            p_initial_condition[2*local_index+1] = 0;
        }


        VecRestoreArray(initial_condition, &p_initial_condition);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);
        //VecView(initial_condition, PETSC_VIEWER_STDOUT_WORLD);
    
        /*
         *  Write data to a file <mOutputFilenamePrefix>_xx.dat, 'xx' refers to 
         *  'xx'th time step using ColumnDataWriter 
         */          
        ///\todo: get writer to write V_m and \phi_e seperately?
        ParallelColumnDataWriter *p_test_writer = NULL;
       
        int time_var_id = 0;
        int voltage_var_id = 0;

        if (mOutputFilenamePrefix.length() > 0)
        {
            mkdir(mOutputDirectory.c_str(), 0777);
                 
            p_test_writer = new ParallelColumnDataWriter(mOutputDirectory,mOutputFilenamePrefix);

            p_test_writer->DefineFixedDimension("Node", "dimensionless", 2*mMesh.GetNumNodes() );
            time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
        
            voltage_var_id = p_test_writer->DefineVariable("Vm_And_Phi_e","mV");
            p_test_writer->EndDefineMode();
        }
        else
        {
            throw Exception("mOutputFilenamePrefix should not be the empty string");
        }

        double current_time = mStartTime;        
        int big_steps = 0; 
        
        p_test_writer->PutVariable(time_var_id, current_time); 
        p_test_writer->PutVector(voltage_var_id, initial_condition);
                        
                        
        while( current_time < mEndTime )
        {
            bidomain_assembler.SetTimes(current_time, current_time+mPdeTimeStep, mPdeTimeStep);
            bidomain_assembler.SetInitialCondition( initial_condition );
            
            mVoltage = bidomain_assembler.Solve(); //(mMesh, mpBidomainPde, bcc);
                                    
            // Free old initial condition
            VecDestroy(initial_condition);

            // Initial condition for next loop is current solution
            initial_condition = mVoltage;
            
            // Writing data out to the file <mOutputFilenamePrefix>.dat                 
            p_test_writer->AdvanceAlongUnlimitedDimension(); //creates a new file
            p_test_writer->PutVariable(time_var_id, current_time); 
            p_test_writer->PutVector(voltage_var_id, mVoltage);

            mpBidomainPde->ResetAsUnsolvedOdeSystem();

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
        //check these have been set
        assert(mLo > -1);
        assert(mHi > -1);        
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
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> & rGetMesh()
    {
        return mMesh;   
    }
    
    BidomainPde<SPACE_DIM>* GetBidomainPde() 
    {
        return mpBidomainPde;  
    }
};


#endif /*BIDOMAINPROBLEM_HPP_*/
