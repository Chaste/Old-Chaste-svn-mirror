#ifndef ABSTRACTCARDIACPROBLEM_HPP_
#define ABSTRACTCARDIACPROBLEM_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "ParallelColumnDataWriter.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "TimeStepper.hpp"
#include "DistributedVector.hpp"
#include "AbstractCardiacPde.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "EventHandler.hpp"
#include "PetscTools.hpp"

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractCardiacProblem
{
private:
    std::string mMeshFilename;
    bool mAllocatedMemoryForMesh;
    
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
    bool mCallChaste2Meshalyzer;
 
    AbstractCardiacPde<SPACE_DIM>* mpCardiacPde;
    
    /** data is not written if output directory or output file prefix are not set*/
    std::string  mOutputDirectory, mOutputFilenamePrefix;

protected:
    AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* mpAssembler; 

    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* mpMesh;
    
    Vec mVoltage; // Current solution
    double mLinearSolverRelativeTolerance;
    double mLinearSolverAbsoluteTolerance;
    bool mUseLinearSolverAbsoluteTolerance;
    

    /**
     * Subclasses must override this method to create a PDE object of the appropriate type.
     * 
     * This class will take responsibility for freeing the object when it is finished with.
     */
    virtual AbstractCardiacPde<SPACE_DIM>* CreateCardiacPde() =0;
    
    /**
     * Subclasses must override this method to create a suitable assembler object.
     * 
     * This class will take responsibility for freeing the object when it is finished with.
     */
    virtual AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* CreateAssembler() =0;

    unsigned mVoltageColumnId;
    unsigned mTimeColumnId;
    unsigned mNodeColumnId;

public: 
    // This (and things in MonodomainProblem) being public are hacks for 
    // AbstractCardiacElectroMechanicsWriter to work.
    // TODO AbstractCardiacElectroMechanicsWriter should be a friend, but not sure
    // how to get friends to work when both friends are templated and abstract.
    ParallelColumnDataWriter* mpWriter;

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
              mpCellFactory(pCellFactory),
              mpMesh(NULL),
              mpWriter(NULL)
    {
        mStartTime        = 0.0;  // ms
        mPdeTimeStep      = 0.01; // ms
        mEndTime          = -1;   // negative so can check has been set
        mPrintingTimeStep = mPdeTimeStep;  // default behaviour: print out every pde time step
        mWriteInfo = false;
        mPrintOutput = true;
        mCallChaste2Meshalyzer = true;
        mpCardiacPde = NULL;
        mpAssembler = NULL;
        mVoltage = NULL;
        mLinearSolverRelativeTolerance=1e-6;
        mUseLinearSolverAbsoluteTolerance = false;
        mAllocatedMemoryForMesh = false;
        
        EventHandler::BeginEvent(EVERYTHING);
    }
    
    virtual ~AbstractCardiacProblem()
    {
        delete mpCardiacPde;
        if (mVoltage)
        {
            VecDestroy(mVoltage);
        }
        
        if(mAllocatedMemoryForMesh)
        {
            delete mpMesh;
        }
        
        EventHandler::EndEvent(EVERYTHING);
    };
    
    /*
     *  Initialise the system. Must be called before Solve()
     */
    void Initialise()
    {
        if (mpMesh==NULL)
        {
            EXCEPTION("SetMesh() or SetMeshFilename() was not set");
        }
        mpCellFactory->SetMesh( mpMesh );
        
        delete mpCardiacPde; // In case we're called twice
        mpCardiacPde = CreateCardiacPde();
    }
    
    void SetLinearSolverRelativeTolerance(const double &rRelTol)
    {
        mLinearSolverRelativeTolerance = rRelTol;
        mUseLinearSolverAbsoluteTolerance = false;        
    }
    
    double GetLinearSolverRelativeTolerance()
    {
        if (mUseLinearSolverAbsoluteTolerance)
        {
            EXCEPTION("No relative tolerance because absolute tolerance set");
        }
        
        return mLinearSolverRelativeTolerance;
    }
    
    void SetLinearSolverAbsoluteTolerance(const double &rAbsTol)
    {
        mLinearSolverAbsoluteTolerance = rAbsTol;
        mUseLinearSolverAbsoluteTolerance = true;        
    }
    
    double GetLinearSolverAbsoluteTolerance()
    {
        if (!mUseLinearSolverAbsoluteTolerance)
        {
            EXCEPTION("No absolute tolerance because relative tolerance set");
        }
        
        return mLinearSolverAbsoluteTolerance;
    }
    
    void PreSolveChecks()
    {
        if ( mpCardiacPde == NULL ) // if pde is NULL, Initialise() probably hasn't been called
        {
            EXCEPTION("Pde is null, Initialise() probably hasn't been called");
        }
        if ( mStartTime > mEndTime )
        {
            EXCEPTION("Start time should be no more than end time");
        }
        if (mPrintOutput==true) 
        {
            if( (mOutputDirectory=="") || (mOutputFilenamePrefix==""))
            {
                EXCEPTION("Either explicitly specify not to print output (call PrintOutput(false)) or specify the output directory and filename prefix");
            }
        }
        if (  mPrintingTimeStep < mPdeTimeStep - 1e-10)
        {
            EXCEPTION("Printing time step is less than Pde time step");
        }
    }
    
    // Perhaps this should be a method of AbstractCardiacPde??)
    
    Vec CreateInitialCondition()
    {
        DistributedVector::SetProblemSize(mpMesh->GetNumNodes());
        Vec initial_condition=DistributedVector::CreateVec(PROBLEM_DIM);
        DistributedVector ic(initial_condition);
        std::vector< DistributedVector::Stripe > stripe;
        stripe.reserve(PROBLEM_DIM);
        
        for (unsigned i=0; i<PROBLEM_DIM; i++)
        {
            stripe.push_back(DistributedVector::Stripe(ic, i));
        }
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            stripe[0][index] = mpCardiacPde->GetCardiacCell(index.Global)->GetVoltage();
            if (PROBLEM_DIM==2)
            {
                stripe[1][index] = 0;
            }
        }
        
        ic.Restore();    
        
        return initial_condition;
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
     *  AFTER setting the pde timestep!
     */
    void PrintEveryNthTimeStep(unsigned n)
    {
        mPrintingTimeStep = n*mPdeTimeStep;
    }
    
    double GetPdeTimeStep()
    {
        return mPdeTimeStep;
    }
    
    
    /** Set whether to call the Chaste2Meshalyzer script.
     * This script gets everything ready to visualize the results with meshalyser
     * and is useful in testing. By default the script is called.
     * In performance testing for example it desirable to disable the script.
     */
    void SetCallChaste2Meshalyzer(bool call)
    {
        mCallChaste2Meshalyzer=call;
    }
    
    void SetMeshFilename(const std::string &rMeshFilename)
    {
        // If this fails the mesh has already been set. We assert rather throw an exception 
        // to avoid a memory leak when checking it throws correctly
        assert(mpMesh==NULL); 
        
        if ( rMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was passed in empty");
        }
        
        mMeshFilename = rMeshFilename;
        
        TrianglesMeshReader<SPACE_DIM, SPACE_DIM> mesh_reader(mMeshFilename);
        mpMesh = new ConformingTetrahedralMesh<SPACE_DIM, SPACE_DIM>();
        mAllocatedMemoryForMesh = true;
        
        EventHandler::BeginEvent(READ_MESH);
        mpMesh->ConstructFromMeshReader(mesh_reader);
        EventHandler::EndEvent(READ_MESH);
    }

    void SetMesh(ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh)
    {
        // If this fails the mesh has already been set. We assert rather throw an exception 
        // to avoid a memory leak when checking it throws correctly
        assert(mpMesh==NULL);
        mAllocatedMemoryForMesh = false;
        assert(pMesh!=NULL);
        mpMesh = pMesh;
    }
    
    void SetOutputDirectory(const std::string& rOutputDirectory)
    {
        mOutputDirectory = rOutputDirectory;
    }
    
    void SetOutputFilenamePrefix(const std::string& rOutputFilenamePrefix)
    {
        mOutputFilenamePrefix = rOutputFilenamePrefix;
    }
    
    /**
     *  Set whether the simulation will generate results files.
     */
    void PrintOutput(bool rPrintOutput)
    {
        mPrintOutput = rPrintOutput;
    }

    /**
     *  Set whether extra info will be written to stdout during computation. 
     */
    void SetWriteInfo(bool writeInfo = true)
    {
        mWriteInfo = writeInfo;
    }
    
    /**
     * Get the final solution vector. This vector is distributed over all processes.
     *
     * In case of Bidomain, this is of length 2*numNodes, and of the form
     *  (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N). 
     *  where V_j is the voltage at node j and phi_j is the
     *  extracellular potential at node j.
     * 
     *  Use with caution since we don't want to alter the state of the PETSc vector.
     */
    Vec GetVoltage()
    {
        return mVoltage;
    }
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> & rGetMesh()    
    {
        return *mpMesh;
    }
    
    AbstractCardiacPde<SPACE_DIM>* GetPde()
    {
        return mpCardiacPde;
    }
    
    void Solve()
    {
        PreSolveChecks();
        mpAssembler = CreateAssembler();
        Vec initial_condition = CreateInitialCondition();
        
//        DistributedVector ic = DistributedVector(initial_condition);
//        DistributedVector::Stripe transmembrane_ic(initial_condition, 0);

        TimeStepper stepper(mStartTime, mEndTime, mPrintingTimeStep);

        if (mPrintOutput)
        {
            EventHandler::BeginEvent(WRITE_OUTPUT);
            InitialiseWriter();
            WriteOneStep(stepper.GetTime(), initial_condition);
            EventHandler::EndEvent(WRITE_OUTPUT);
        }
        
        // If we have already run a simulation, free the old solution vec
        if (mVoltage)
        {
            VecDestroy(mVoltage);
        }
        
        while ( !stepper.IsTimeAtEnd() )
        {
            // solve from now up to the next printing time
            mpAssembler->SetTimes(stepper.GetTime(), stepper.GetNextTime(), mPdeTimeStep);
            mpAssembler->SetInitialCondition( initial_condition );

            try
            {
                mVoltage = mpAssembler->Solve();
            }
            catch (Exception &e)
            {
                // Free memory.
                delete mpAssembler;                 
                //VecDestroy(initial_condition);
                // Close files
                if (mPrintOutput)
                {
                    mpWriter->Close();
                    delete mpWriter;
                }
                PetscTools::ReplicateException(true);
                // Re-throw
                throw e;
            }
            PetscTools::ReplicateException(false);
            
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
                EventHandler::BeginEvent(WRITE_OUTPUT);
                mpWriter->AdvanceAlongUnlimitedDimension(); //creates a new file
                WriteOneStep(stepper.GetTime(), mVoltage);
                EventHandler::EndEvent(WRITE_OUTPUT);
            }
        }
        
        // Free assembler
        delete mpAssembler;

        // close the file that stores voltage values
        if (mPrintOutput)
        {
            bool am_master = mpWriter->AmMaster();
            mpWriter->Close();
            delete mpWriter;
            
            if (am_master && mCallChaste2Meshalyzer) // ie only if master process and results files were written
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
    
    virtual void DefineWriterColumns()
    {
        mNodeColumnId = mpWriter->DefineFixedDimension("Node", "dimensionless", mpMesh->GetNumNodes() );
        mTimeColumnId = mpWriter->DefineUnlimitedDimension("Time","msecs");
        mVoltageColumnId = mpWriter->DefineVariable("V","mV");
    }
    
    virtual void WriteOneStep(double time, Vec voltageVec)
    {
        if (mpWriter->AmMaster())
        {
            for (unsigned node = 0; node < mpMesh->GetNumNodes(); node++)
            {
                mpWriter->PutVariable(mNodeColumnId, node, node);
            }
            mpWriter->PutVariable(mTimeColumnId, time);
        }
        DistributedVector::Stripe transmembrane(voltageVec, 0);
        mpWriter->PutVectorStripe(mVoltageColumnId, transmembrane);
    }
    
    void InitialiseWriter()
    {
        mpWriter = new ParallelColumnDataWriter(mOutputDirectory,mOutputFilenamePrefix);
        DefineWriterColumns();
        mpWriter->EndDefineMode();
    }        
};

#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
