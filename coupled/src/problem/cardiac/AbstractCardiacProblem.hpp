#ifndef ABSTRACTCARDIACPROBLEM_HPP_
#define ABSTRACTCARDIACPROBLEM_HPP_



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
       
    }
};

#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
