#ifndef MONODOMAINPROBLEM_HPP_
#define MONODOMAINPROBLEM_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.cpp"
#include "ParallelColumnDataWriter.hpp"
#include "MonodomainPde.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"

/**
 * Class which specifies and solves a monodomain problem.
 */
template<unsigned SPACE_DIM>
class MonodomainProblem
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
    
    MonodomainPde<SPACE_DIM> *mpMonodomainPde;
    
    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    
    Vec mVoltage; // Current solution
    unsigned mLo, mHi;
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> mMesh;
    
public:

    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should 
     * create cells.
     */
    MonodomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            : mMeshFilename(""),      // i.e. undefined
            mOutputDirectory(""),   // i.e. undefined
            mOutputFilenamePrefix(""),   // i.e. undefined
            mpMonodomainPde(NULL),
            mpCellFactory(pCellFactory)
    {
        mStartTime        = 0.0;  // ms
        mPdeTimeStep      = 0.01; // ms
        mEndTime          = -1;   // negative so can check has been set
        mPrintingTimeStep = -1;   // negative so can check has been set
        
        mWriteInfo = false;
    }
    
    /**
     * Destructor
     */
    ~MonodomainProblem()
    {
        delete mpMonodomainPde;
    }
    
    /** Initialise the system. Must be called before Solve() */
    void Initialise()
    {
        if ( mMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was passed in empty");
        }
        
        TrianglesMeshReader<SPACE_DIM, SPACE_DIM> mesh_reader(mMeshFilename);
        mMesh.ConstructFromMeshReader(mesh_reader);
        
        mpCellFactory->SetMesh( &mMesh );
        
        mpMonodomainPde = new MonodomainPde<SPACE_DIM>(mpCellFactory);
    }
    
    /**
     * Solve the problem
     */
    void Solve()
    {
        if ( mpMonodomainPde == NULL ) // if pde is NULL, Initialise() probably hasn't been called
        {
            EXCEPTION("Monodomain pde is null, Initialise() probably hasn't been called");
        }
        
        if ( mStartTime >= mEndTime )
        {
            EXCEPTION("Start time should be less than end time");
        }
        
        // Assembler
        MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM> monodomain_assembler(&mMesh, mpMonodomainPde);
        
        // initial condition;
        Vec initial_condition= DistributedVector::CreateVec();
        DistributedVector ic(initial_condition);
        mLo = DistributedVector::Begin().Global;
        mHi = DistributedVector::End().Global;
        // Set a constant initial voltage throughout the mMesh
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
             {
                ic[index] = mpMonodomainPde->GetCardiacCell(index.Global)->GetVoltage();
             }
        ic.Restore();        


        
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
            
            p_test_writer->DefineFixedDimension("Node", "dimensionless", mMesh.GetNumNodes() );
            time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
            
            voltage_var_id = p_test_writer->DefineVariable("V","mV");
            p_test_writer->EndDefineMode();
        }
        
        double current_time = mStartTime;
        
        if (write_files)
        {
            p_test_writer->PutVariable(time_var_id, current_time);
            p_test_writer->PutVector(voltage_var_id, initial_condition);
        }
        
        
        if ( mPrintingTimeStep < 0) //ie if it hasn't been set
        {
            mPrintingTimeStep = mPdeTimeStep; ///\todo: pick good default
        }
        
        
// at the moment the tend-tstart must be a multiple of the printing timestep
// because the 'next_printing_time = mEndTime' line below has been commented
// out (see comment below). remove this assert when this is sorted out
        assert( fabs(        ( (mEndTime-mStartTime)/mPrintingTimeStep )
                             - round( (mEndTime-mStartTime)/mPrintingTimeStep ) ) < 1e-10 );
                             
                             
        // check the printing time step is a multiple of the pde timestep.
        assert(  fabs( (mPrintingTimeStep/mPdeTimeStep)
                       -round(mPrintingTimeStep/mPdeTimeStep) ) < 1e-10 );
                       
        while ( current_time < mEndTime )
        {
            // compute the next printing time
            double next_printing_time = current_time + mPrintingTimeStep;
            if (next_printing_time > mEndTime)
            {
// this line is needed but it's use leads to assertions tripping due to
// floating point errors, so it's been commented out. see ticket 152
///\todo: sort this out! then change TestPrintsOnlyAtRequestedTimes in
// TestMonodomainDg0Assembler and TestBidomainProblem to have an endtime that
// is not a multiple of the printing time...
//                next_printing_time = mEndTime;
            }
            
            // solve from now up to the next printing time
            monodomain_assembler.SetTimes(current_time, next_printing_time, mPdeTimeStep);
            monodomain_assembler.SetInitialCondition( initial_condition );
            
            try
            {
                mVoltage = monodomain_assembler.Solve();
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
            if (mWriteInfo)
            {
                WriteInfo(current_time);
            }
            
            // Writing data out to the file <mOutputFilenamePrefix>.dat
            if (write_files)
            {
                p_test_writer->AdvanceAlongUnlimitedDimension(); // creates a new file
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
        
        PetscInt my_rank;
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
        if (pdeTimeStep <= 0)
        {
            EXCEPTION("Pde time step should be positive");
        }
        mPdeTimeStep = pdeTimeStep;
    }
    
    /** Set the times to print output. The printing time step must be
     *  a multiple of the pde timestep 
     */
    void SetPrintingTimeStep(double printingTimeStep)
    {
        if (printingTimeStep <= 0.0)
        {
            EXCEPTION("Printing time step should be positive");
        }
        mPrintingTimeStep = printingTimeStep;
    }
    
    /** Set the simulation to print every n timesteps. Only set this
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
    
    /**
     * Get the final solution vector. This vector is distributed over all processes.
     */
    
    Vec GetVoltage()
    {
        //Use with caution since we don't want to alter the state of the PETSc vector
        return mVoltage;
    }
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> & rGetMesh()
    {
        return mMesh;
    }
    
    MonodomainPde<SPACE_DIM> * GetMonodomainPde()
    {
        return mpMonodomainPde;
    }
    
    /**
     *  Set info to be printed during computation. 
     */
    void SetWriteInfo(bool writeInfo = true)
    {
        mWriteInfo = writeInfo;
    }
    
    
    /**
     *  Print out time and max/min voltage values at current time.
     
     */
    void WriteInfo(double time)
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
        
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(mVoltage);
        
        double v_max = -1e5, v_min = 1e5;
        for (unsigned i=0; i<mMesh.GetNumNodes(); i++)
        {
            if ( voltage_replicated[i] > v_max)
            {
                v_max = voltage_replicated[i];
            }
            if ( voltage_replicated[i] < v_min)
            {
                v_min = voltage_replicated[i];
            }
        }
        
        std::cout << " max/min V = "
        <<   v_max << " "
        <<   v_min << "\n" << std::flush;
    }
};

#endif /*MONODOMAINPROBLEM_HPP_*/
