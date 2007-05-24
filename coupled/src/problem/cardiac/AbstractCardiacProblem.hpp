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
    
    Vec CreateInitialCondition(AbstractCardiacPde<SPACE_DIM>* pCardiacPde, unsigned num_domains)
    {
        assert(num_domains<=2);
        DistributedVector::SetProblemSize(this->mMesh.GetNumNodes());
        Vec initial_condition=DistributedVector::CreateVec(num_domains);
        DistributedVector ic(initial_condition);
        
        std::vector< DistributedVector::Stripe > stripe;
        
        for (unsigned i=0; i<num_domains; i++)
        {
            stripe.push_back(DistributedVector::Stripe(ic, i));
        }
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            stripe[0][index]= pCardiacPde->GetCardiacCell(index.Global)->GetVoltage();
            
            if (num_domains==2)
            {
                stripe[1][index] =0;
            }
        }
        
        ic.Restore();
        
        return initial_condition;
    }
    
};

#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
