#ifndef ABSTRACTCARDIACPROBLEM_HPP_
#define ABSTRACTCARDIACPROBLEM_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "ParallelColumnDataWriter.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "TimeStepper.hpp"
#include "DistributedVector.hpp"

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractCardiacProblem
{
protected:
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
    bool mPrintOutput;
    
    /** data is not written if output directory or output file prefix are not set*/
    std::string  mOutputDirectory, mOutputFilenamePrefix;
        
    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    Vec mVoltage; // Current solution
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> mMesh;

public:    
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should 
     * create cells.
     */
    AbstractCardiacProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            : mMeshFilename(""),     // i.e. undefined
            mOutputDirectory(""),  // i.e. undefined
            mOutputFilenamePrefix(""),   // i.e. undefined
            mpCellFactory(pCellFactory)
    {
        mStartTime        = 0.0;  // ms
        mPdeTimeStep      = 0.01; // ms
        mEndTime          = -1;   // negative so can check has been set
        mPrintingTimeStep = mPdeTimeStep;  // default behaviour: print out every pde time step
        mWriteInfo = false;
        mPrintOutput = true;
    }
    
    virtual ~AbstractCardiacProblem() {};
    
    /*
     *  Initialise the system. Must be called before Solve()
     */
    void Initialise(AbstractCardiacPde<SPACE_DIM>* pCardiacPde)
    {
        if ( mMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was not set");
        }
        mpCellFactory->SetMesh( &mMesh );
        if (pCardiacPde)
        {
            delete pCardiacPde;
        }
    }
    
    void PreSolveChecks(AbstractCardiacPde<SPACE_DIM>* pCardiacPde)
    {
        if ( pCardiacPde == NULL ) // if pde is NULL, Initialise() probably hasn't been called
        {
            EXCEPTION("Pde is null, Initialise() probably hasn't been called");
        }
        if ( mStartTime >= mEndTime )
        {
            EXCEPTION("Start time should be less than end time");
        }
    }
    
    // Perhaps this should be a method of AbstractCardiacPde??)
    
    Vec CreateInitialCondition(AbstractCardiacPde<SPACE_DIM>* pCardiacPde)
    {
        DistributedVector::SetProblemSize(mMesh.GetNumNodes());
        Vec initial_condition=DistributedVector::CreateVec(PROBLEM_DIM);
        DistributedVector ic(initial_condition);
        std::vector< DistributedVector::Stripe > stripe;
        
        for (unsigned i=0; i<PROBLEM_DIM; i++)
        {
            stripe.push_back(DistributedVector::Stripe(ic, i));
        }
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            stripe[0][index]= pCardiacPde->GetCardiacCell(index.Global)->GetVoltage();
            if (PROBLEM_DIM==2)
            {
                stripe[1][index] =0;
            }
        }
        
        ic.Restore();    
        
        return initial_condition;
    }
    
    /**
     *  Set the simulation to print every n timesteps. Only set this
     *  AFTER setting the pde timestep
     */
    void PrintOutput(const bool& rPrintOutput)
    {
        mPrintOutput = rPrintOutput;
    }
    
    void SetStartTime(const double &rStartTime)
    {
        mStartTime = rStartTime;
    }
    
    void SetEndTime(const double &rEndTime)
    {
        mEndTime = rEndTime;
    }
    
    /**
     * Set the PDE time step.
     * \todo SetPdeAndPrintingTimeStep
     * Note that the printing time step should also set with this call.
     */
    void SetPdeTimeStep(double pdeTimeStep)
    {
        if (pdeTimeStep <= 0)
        {
            EXCEPTION("Pde time step should be positive");
        }
        mPdeTimeStep = pdeTimeStep;
    }
    
    /** 
     * Set the times to print output. The printing time step must be
     * a multiple of the pde timestep. If SetPdeTimeStep is used it should be called
     * before SetPrintingTimeStep.
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
        if ( mMeshFilename!="" )
        {
            EXCEPTION("Mesh filename was already set");
        }
        if ( rMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was passed in empty");
        }
        mMeshFilename = rMeshFilename;
        
        TrianglesMeshReader<SPACE_DIM, SPACE_DIM> mesh_reader(mMeshFilename);
        mMesh.ConstructFromMeshReader(mesh_reader);
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
     *
     * In case of Bidomain, this is of length 2*numNodes, and of the form
     *  (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N). 
     *  where V_j is the voltage at node j and phi_j is the
     *  extracellular potential at node j.
     * 
     *  This vector is distributed over all processes.
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

    /**
     *  Set info to be printed during computation. 
     */
    void SetWriteInfo(bool writeInfo = true)
    {
        mWriteInfo = writeInfo;
    }
    
    void Solve(AbstractLinearDynamicProblemAssembler<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>& assembler,
               AbstractCardiacPde<SPACE_DIM>* pCardiacPde)
    {
        Vec initial_condition = AbstractCardiacProblem<SPACE_DIM, PROBLEM_DIM>::CreateInitialCondition(pCardiacPde);
        ParallelColumnDataWriter *p_test_writer = NULL;
        unsigned time_var_id = 0;
        unsigned voltage_var_id = 0;
        bool write_files = false;

        TimeStepper stepper(mStartTime, mEndTime, mPrintingTimeStep);

        if (mPrintOutput)
        {
            write_files = true;
            p_test_writer = new ParallelColumnDataWriter(mOutputDirectory,mOutputFilenamePrefix);
            p_test_writer->DefineFixedDimension("Node", "dimensionless", PROBLEM_DIM*mMesh.GetNumNodes() );
            time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
            voltage_var_id = p_test_writer->DefineVariable(ColumnName(),"mV");
            p_test_writer->EndDefineMode();
            p_test_writer->PutVariable(time_var_id, stepper.GetTime());
            p_test_writer->PutVector(voltage_var_id, initial_condition);
        }
        
        while ( !stepper.IsTimeAtEnd() )
        {
            // solve from now up to the next printing time
            assembler.SetTimes(stepper.GetTime(), stepper.GetNextTime(), mPdeTimeStep);
            assembler.SetInitialCondition( initial_condition );

            try
            {
                mVoltage = assembler.Solve();
            }
            //Ill-conditioned solutions are covered in Monodomain problem
            //(and possibly in Nightly/Weekly) so we don't insist on it
            //in the coverage test.
            #define COVERAGE_IGNORE
            catch (Exception &e)
            {
                if (mPrintOutput)
                {
                    p_test_writer->Close();
                    delete p_test_writer;
                }
                throw e;
            }
            #undef COVERAGE_IGNORE
            
            // Free old initial condition
            VecDestroy(initial_condition);
            
            // Initial condition for next loop is current solution
            initial_condition = mVoltage;
            
            // update the current time
            stepper.AdvanceOneTimeStep();
            
            if (mPrintOutput)
            {
                // print out details at current time if asked for
                if (mWriteInfo)
                {
                    WriteInfo(stepper.GetTime());
                }
                
                // Writing data out to the file <mOutputFilenamePrefix>.dat
                p_test_writer->AdvanceAlongUnlimitedDimension(); //creates a new file
                p_test_writer->PutVariable(time_var_id, stepper.GetTime());
                p_test_writer->PutVector(voltage_var_id, mVoltage);
            }
        }

        // close the file that stores voltage values
        if (mPrintOutput)
        {

            p_test_writer->Close();
            delete p_test_writer;

            
            PetscInt my_rank;
            MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
            if (my_rank==0) // ie only if master process and results files were written
            {
                // call shell script which converts the data to meshalyzer format
                std::string chaste_2_meshalyzer;
                std::stringstream space_dim;
                space_dim << SPACE_DIM;
                chaste_2_meshalyzer = "anim/chaste2meshalyzer "     // the executable.
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
    }
    
    virtual void WriteInfo(double time) =0;
    
    virtual std::string ColumnName() =0;
};

#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
