#ifndef ABSTRACTCARDIACPROBLEM_HPP_
#define ABSTRACTCARDIACPROBLEM_HPP_

#include "DistributedVector.hpp"


template<unsigned SPACE_DIM>
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
    
    // perhaps should be stored in the cardiac pde
    unsigned mNumDomains;
    
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
    
    void Initialise(AbstractCardiacPde<SPACE_DIM>* pCardiacPde)
    {
        if ( this->mMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was not set");
        }
        
        this->mpCellFactory->SetMesh( &this->mMesh );
        
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
        
        if ( this->mStartTime >= this->mEndTime )
        {
            EXCEPTION("Start time should be less than end time");
        }
    }
    
    // Perhaps this should be a method of AbstractCardiacPde??)
    
    Vec CreateInitialCondition(AbstractCardiacPde<SPACE_DIM>* pCardiacPde)
    {
        DistributedVector::SetProblemSize(this->mMesh.GetNumNodes());
        Vec initial_condition=DistributedVector::CreateVec(mNumDomains);
        DistributedVector ic(initial_condition);
        
        std::vector< DistributedVector::Stripe > stripe;
        
        for (unsigned i=0; i<mNumDomains; i++)
        {
            stripe.push_back(DistributedVector::Stripe(ic, i));
        }
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            stripe[0][index]= pCardiacPde->GetCardiacCell(index.Global)->GetVoltage();
            
            if (mNumDomains==2)
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
        this->mPrintOutput = rPrintOutput;
    }
    
        void SetStartTime(const double &rStartTime)
    {
        this->mStartTime = rStartTime;
    }
    
    void SetEndTime(const double &rEndTime)
    {
        this->mEndTime = rEndTime;
    }
    
    /**
     * Set the PDE time step.
     * \todo SetPdeAndPrintingTimeStep
     * Note that the printing time step should also set with this call.
     * Move assertion
     */
    
    void SetPdeTimeStep(double pdeTimeStep)
    {
        if (pdeTimeStep <= 0)
        {
            EXCEPTION("Pde time step should be positive");
        }
        this->mPdeTimeStep = pdeTimeStep;
    }
    
    /** Set the times to print output. The printing time step must be
     *  a multiple of the pde timestep. If SetPdeTimeStep is used it should be called
     * before SetPrintingTimeStep.
     */
    void SetPrintingTimeStep(double printingTimeStep)
    {
        if (printingTimeStep <= 0.0)
        {
            EXCEPTION("Printing time step should be positive");
        }
        this->mPrintingTimeStep = printingTimeStep;
    }
    
    /** Set the simulation to print every n timesteps. Only set this
     *  AFTER setting the pde timestep
     */
    void PrintEveryNthTimeStep(unsigned n)
    {
        this->mPrintingTimeStep = n*this->mPdeTimeStep;
    }
    
    
    double GetPdeTimeStep()
    {
        return this->mPdeTimeStep;
    }
    
    void SetMeshFilename(const std::string &rMeshFilename)
    {
        if ( this->mMeshFilename!="" )
        {
            EXCEPTION("Mesh filename was already set");
        }
        if ( rMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was passed in empty");
        }
        this->mMeshFilename = rMeshFilename;
        
        TrianglesMeshReader<SPACE_DIM, SPACE_DIM> mesh_reader(this->mMeshFilename);
        this->mMesh.ConstructFromMeshReader(mesh_reader);
    }
    
    void SetOutputDirectory(const std::string &rOutputDirectory)
    {
        this->mOutputDirectory = rOutputDirectory;
    }
    
    void SetOutputFilenamePrefix(const std::string &rOutputFilenamePrefix)
    {
        this->mOutputFilenamePrefix = rOutputFilenamePrefix;
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
        return this->mVoltage;
    }
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> & rGetMesh()
    {
        return this->mMesh;
    }

    
    /**
     *  Set info to be printed during computation. 
     */
    void SetWriteInfo(bool writeInfo = true)
    {
        this->mWriteInfo = writeInfo;
    }
    
    
};

#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
