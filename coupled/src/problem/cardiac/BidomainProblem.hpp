#ifndef BIDOMAINPROBLEM_HPP_
#define BIDOMAINPROBLEM_HPP_


#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "BoundaryConditionsContainer.hpp"
#include "BidomainDg0Assembler.hpp"
#include "TrianglesMeshReader.cpp"
#include "ParallelColumnDataWriter.hpp"
#include "BidomainPde.hpp"
#include "AbstractCardiacCellFactory.hpp"


/**
 * Class which specifies and solves a bidomain problem.
 * 
 * The solution vector is of the form:
 * (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N), 
 * where V_j is the voltage at node j and phi_j is the 
 * extracellular potential at node j.
 */
template<int SPACE_DIM>
class BidomainProblem
{
private:
    std::string mMeshFilename;

    /** 
     *  Start time defaults to 0, pde timestep defaults to 0.01 (ms), the
     *  end time is not defaulted and must be set
     */
    double mStartTime;
    double mEndTime;
    double mPdeTimeStep;
    double mPrintingTimeStep;       
    
    bool mWriteInfo; 

    /** data is not written if output directory or output file prefix are not set*/ 
    std::string  mOutputDirectory, mOutputFilenamePrefix;  

    BidomainPde<SPACE_DIM>* mpBidomainPde;

    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    
    Vec mVoltage; // Current solution
    unsigned mLo, mHi; // Ownership range for the double size potential vector mVoltage

    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> mMesh;
    
public:
    
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should 
     * create cells.
     */
    BidomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
    : mMeshFilename(""),     // i.e. undefined
      mOutputDirectory(""),  // i.e. undefined
      mOutputFilenamePrefix(""),   // i.e. undefined
      mpBidomainPde(NULL),
      mpCellFactory(pCellFactory)
    {
        mStartTime        = 0.0;  // ms
        mPdeTimeStep      = 0.01; // ms
        mEndTime          = -1;   // negative so can check has been set
        mPrintingTimeStep = -1;   // negative so can check has been set

        //initialise these to 0
        mLo = 0;
        mHi = 0; 
        
        mWriteInfo = false;
    }

    /**
     * Destructor
     */
    ~BidomainProblem()
    { 
        delete mpBidomainPde;
    }
    
    
    /** Initialise the system. Must be called before Solve() */
    void Initialise()
    {
        assert( mMeshFilename!="" );

        TrianglesMeshReader<SPACE_DIM, SPACE_DIM> mesh_reader(mMeshFilename);
        mMesh.ConstructFromMeshReader(mesh_reader);
        
        mpCellFactory->SetMesh( &mMesh );
        
        mpBidomainPde = new BidomainPde<SPACE_DIM>( mpCellFactory, mPdeTimeStep );
    }
     

    /**
     * Solve the problem
     */
    void Solve()
    {
        assert( mpBidomainPde != NULL );  // if pde is NULL, Initialise() probably hasn't been called
        assert( mStartTime < mEndTime );
        
        // Linear solver
        SimpleLinearSolver linear_solver;
    
        // Assembler
        BidomainDg0Assembler<SPACE_DIM,SPACE_DIM> bidomain_assembler(mpBidomainPde, &mMesh, &linear_solver);
        
        // initial condition;   
        Vec initial_condition;
        
        unsigned big_lo, big_hi;
        mpBidomainPde->GetOwnershipRange(big_lo, big_hi);

        VecCreateMPI(PETSC_COMM_WORLD, 2*(big_hi-big_lo) , 2*mMesh.GetNumNodes(), &initial_condition);  
        double* p_initial_condition;
        VecGetArray(initial_condition, &p_initial_condition); 

        PetscInt temp_lo, temp_hi;
        VecGetOwnershipRange(initial_condition, &temp_lo, &temp_hi);
        mLo=(unsigned) temp_lo;
        mHi=(unsigned) temp_hi;
        
        for (unsigned global_index=big_lo; global_index<big_hi; global_index++)
        {
            unsigned local_index = global_index - big_lo;
            p_initial_condition[2*local_index  ] = mpBidomainPde->GetCardiacCell(global_index)->GetVoltage();
            p_initial_condition[2*local_index+1] = 0;
        }


        VecRestoreArray(initial_condition, &p_initial_condition);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);
        //VecView(initial_condition, PETSC_VIEWER_STDOUT_WORLD);
    
        //  Write data to a file <mOutputFilenamePrefix>_xx.dat, 'xx' refers to 
        //  'xx'th time step using ColumnDataWriter         
        ParallelColumnDataWriter *p_test_writer = NULL;
       
        unsigned time_var_id = 0;
        unsigned voltage_var_id = 0;
        bool write_files = false;

        if (mOutputFilenamePrefix.length() > 0)
        {
            write_files = true;

            p_test_writer = new ParallelColumnDataWriter(mOutputDirectory,mOutputFilenamePrefix);

            p_test_writer->DefineFixedDimension("Node", "dimensionless", 2*mMesh.GetNumNodes() );
            time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
        
            voltage_var_id = p_test_writer->DefineVariable("Vm_And_Phi_e","mV");
            p_test_writer->EndDefineMode();
        }

        double current_time = mStartTime;        

        if (write_files)
        {        
            p_test_writer->PutVariable(time_var_id, current_time); 
            p_test_writer->PutVector(voltage_var_id, initial_condition);
        }     
 
        if( mPrintingTimeStep < 0) //ie if it hasn't been set
        {
            mPrintingTimeStep = mPdeTimeStep; ///\todo: pick good default
        }    
        
        
// at the moment the tend-tstart must be a multiple of the printing timestep 
// because the 'next_printing_time = mEndTime' line below has been commented
// out (see comment below). remove this assert when this is sorted out
        assert( fabs(        ( (mEndTime-mStartTime)/mPrintingTimeStep )
                      - round( (mEndTime-mStartTime)/mPrintingTimeStep ) ) < 1e-10 );     
        
        
        // check the printing time step is a multiple of the pde timestep.
        assert( fabs(        (mPrintingTimeStep/mPdeTimeStep)
                       -round(mPrintingTimeStep/mPdeTimeStep) ) < 1e-10 );   
 
                        
        while( current_time < mEndTime )
        {
            // compute the next printing time
            double next_printing_time = current_time + mPrintingTimeStep;
            if(next_printing_time > mEndTime)
            {
// this line is needed but it's use leads to assertions tripping due to
// floating point errors. 
///\todo: sort this out! then change TestPrintsOnlyAtRequestedTimes in 
// TestMonodomainDg0Assembler and TestBidomainProblem to have an endtime that
// is not a multiple of the printing time...
//                next_printing_time = mEndTime;
            }
            
            // solve from now up to the next printing time
            bidomain_assembler.SetTimes(current_time, next_printing_time, mPdeTimeStep);
            bidomain_assembler.SetInitialCondition( initial_condition );
           
            try
            { 
                mVoltage = bidomain_assembler.Solve();
            } 
            catch (Exception &e) 
            {
                if (write_files)
                {
                    p_test_writer->Close();
                    delete p_test_writer;
                }
                throw e;
            }
                                    
            // Free old initial condition
            VecDestroy(initial_condition);

            // Initial condition for next loop is current solution
            initial_condition = mVoltage;
            
            // update the current time
            current_time = next_printing_time; 
            
            // print out details at current time if asked for
            if(mWriteInfo)
            {
                WriteInfo(current_time);
            } 
            
            // Writing data out to the file <mOutputFilenamePrefix>.dat                 
            if (write_files)
            {
                p_test_writer->AdvanceAlongUnlimitedDimension(); //creates a new file
                p_test_writer->PutVariable(time_var_id, current_time); 
                p_test_writer->PutVector(voltage_var_id, mVoltage);
            }
        }

        // close the file that stores voltage values            
        if (write_files)
        {
            p_test_writer->Close();
            delete p_test_writer;
        }
        
        int my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if ((my_rank==0) && (write_files)) // ie only if master process and results files were written
        {
            // call shell script which converts the data to meshalyzer format
            std::string chaste_2_meshalyzer;
            std::stringstream space_dim;
            space_dim << SPACE_DIM;
            chaste_2_meshalyzer = "anim/chaste2meshalyzer "         // the executable.
                                      + space_dim.str() + " "       // argument 1 is the dimension.
                                      + mMeshFilename + " "         // arg 2 is mesh prefix, path relative to 
                                                                    // the main chaste directory.
                                      + mOutputDirectory + "/" 
                                      + mOutputFilenamePrefix + " " // arg 3 is the results folder and prefix, 
                                                                    // relative to the testoutput folder. 
                                      + "last_simulation";          // arg 4 is the output prefix, relative to 
                                                                    // anim folder.
                                            
            system(chaste_2_meshalyzer.c_str());
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
        mPdeTimeStep = pdeTimeStep;
    }
    
    /** 
     *  Set the times to print output. The printing time step must be 
     *  a multiple of the pde timestep 
     */
    void SetPrintingTimeStep(double printingTimeStep)
    {
        assert(0.0 < printingTimeStep);
        mPrintingTimeStep = printingTimeStep;
    }
    
    /** 
     *  Set the simulation to print every n timesteps. Only set this
     *  AFTER setting the pde timestep
     */
    void PrintEveryNthTimeStep(unsigned n)
    {
        mPrintingTimeStep = n*mPdeTimeStep;
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
    
 
    /** Get the final solution vector. This is of length 2*numNodes, and of the form
     *  (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N). 
     *  where V_j is the voltage at node j and phi_j is the
     *  extracellular potential at node j.
     * 
     *  This vector is distributed over all processes,
     *  with the current process owning the [lo, ..., hi-1] components of the vector.
     */
    void GetVoltageArray(double **pVoltageArray, unsigned &lo, unsigned &hi)
    {
        //check these have been set
        assert(mHi > 0);        
        VecGetArray(mVoltage, pVoltageArray);
        lo=mLo;
        hi=mHi; 
    }
    
    Vec GetVoltage()
    {
        //Use with caution since we don't want to alter the state of the PETSc vector
        return mVoltage;
    }
    
    /** call this after GetVoltageArray to avoid memory leaks*/
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
    
    /** 
     *  Get the pde. Can only be called after Initialise()
     */
    BidomainPde<SPACE_DIM>* GetBidomainPde() 
    {
        assert(mpBidomainPde!=NULL);
        return mpBidomainPde;  
    }
    
    /** 
     *  Set info to be printed during computation. 
     */
    void SetWriteInfo(bool writeInfo = true)
    {
        mWriteInfo = writeInfo;
    }
    
    
    /**
     *  Print out time and max/min voltage/phi_e values at current time.
     */
    void WriteInfo(double time)
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
                 
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(mVoltage);
        
        double v_max = -1e5, v_min = 1e5, phi_max = -1e5, phi_min = 1e5;
        for(int i=0; i<mMesh.GetNumNodes(); i++)
        {
            if( voltage_replicated[2*i] > v_max)
            {
                v_max = voltage_replicated[2*i];
            }
            if( voltage_replicated[2*i] < v_min)
            {
               v_min = voltage_replicated[2*i];
            }
            if( voltage_replicated[2*i+1] > phi_max)
            {
                phi_max = voltage_replicated[2*i+1];
            }
            if( voltage_replicated[2*i+1] < phi_min)                    
            {
                phi_min = voltage_replicated[2*i+1];
            }

        }
        
        std::cout << " max/min V, phi_e = " 
                  << v_max << " "  
                  << v_min << " "  
                  << phi_max << " "  
                  << phi_min << "\n" << std::flush;  
    }
};


#endif /*BIDOMAINPROBLEM_HPP_*/
