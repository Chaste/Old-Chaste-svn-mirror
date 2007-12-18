#include "TissueCell.hpp"
#include "CellTypes.hpp"
#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"


TissueCell::TissueCell(CellType cellType,
                       CellMutationState mutationState,
                       AbstractCellCycleModel *pCellCycleModel,
                       bool archiving)
    : mpCellCycleModel(pCellCycleModel)
{
    if (SimulationTime::Instance()->IsStartTimeSetUp()==false)
    {
        EXCEPTION("TissueCell is setting up a cell cycle model but SimulationTime has not been set up");
    }
    if(pCellCycleModel==NULL)
    {
        EXCEPTION("Cell cycle model is null");
    }
    
    mCellType = cellType;
    mMutationState = mutationState;
    mCanDivide = false;
    mUndergoingApoptosis = false;
    mIsDead = false;
    mDeathTime = DBL_MAX; // this has to be initialised for archiving...
    mNodeIndex = (unsigned)(-1); // initialise to a silly value for archiving (avoid memory check error)
    mIsLogged = false;
    if (archiving)
    {	// If we are called by the archiver ONLY set the cell (as in abstract class)
    	// don't reset the cell cycle variables as the standard set cell might.
    	mpCellCycleModel->AbstractCellCycleModel::SetCell(this);
    }
    else
    {
    	mpCellCycleModel->SetCell(this);
    }
    mSymmetricDivision = false;
    mAncestor = 0u; // Has to be set by a SetAncestor() call (usually from Tissue)
}

void TissueCell::CommonCopy(const TissueCell &other_cell)
{
    // Copy private data members
    mCanDivide = other_cell.mCanDivide;
    // Copy 'easy' protected data members
    mCellType = other_cell.mCellType;
    mMutationState = other_cell.mMutationState;
    mUndergoingApoptosis = other_cell.mUndergoingApoptosis;
    mIsDead = other_cell.mIsDead;
    mDeathTime = other_cell.mDeathTime;
    mNodeIndex = other_cell.mNodeIndex;
    mIsLogged = other_cell.mIsLogged;
    mSymmetricDivision = other_cell.mSymmetricDivision;
    mAncestor = other_cell.mAncestor;
    
    // Copy cell cycle model
    // First create a new object
    mpCellCycleModel = other_cell.mpCellCycleModel->CreateCellCycleModel();
    // Then copy its state
    *mpCellCycleModel = *(other_cell.mpCellCycleModel);
    // note: we call the base class version because we want to do model.mpCell=*this
    // only, as the model is fully set up (from the above line) already.
    mpCellCycleModel->AbstractCellCycleModel::SetCell(this);
}


TissueCell::TissueCell(const TissueCell &other_cell)
{
    CommonCopy(other_cell);
}


TissueCell& TissueCell::operator=(const TissueCell &other_cell)
{
    // In case this is self-assignment, don't delete the cell cycle model...
    AbstractCellCycleModel* temp = mpCellCycleModel;
    CommonCopy(other_cell);
    // ...until after we've copied it.
    delete temp;
    return *this;
}


TissueCell::~TissueCell()
{
    delete mpCellCycleModel;
}


void TissueCell::SetBirthTime(double birthTime)
{
    mpCellCycleModel->SetBirthTime(birthTime);
}


void TissueCell::SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel)
{
    if (mpCellCycleModel != pCellCycleModel)
    {
        delete mpCellCycleModel;
    }
    mpCellCycleModel = pCellCycleModel;
    mpCellCycleModel->SetCell(this);
}


AbstractCellCycleModel* TissueCell::GetCellCycleModel() const
{
    return mpCellCycleModel;
}


void TissueCell::InitialiseCellCycleModel()
{
    mpCellCycleModel->Initialise();
}


/**
 * Set the node at which this cell is positioned.
 * 
 * @param index Index of the node
 */
void TissueCell::SetNodeIndex(unsigned index)
{
    mNodeIndex = index;
}


unsigned TissueCell::GetNodeIndex() const
{
    return mNodeIndex;
}


double TissueCell::GetAge() const
{
    return mpCellCycleModel->GetAge();
}


double TissueCell::GetBirthTime() const
{
    return mpCellCycleModel->GetBirthTime();
}


CellType TissueCell::GetCellType() const
{
    return mCellType;
}


CellMutationState TissueCell::GetMutationState() const
{
    return mMutationState;
}


void TissueCell::SetCellType(CellType cellType)
{
    mCellType = cellType;
}


void TissueCell::SetMutationState(CellMutationState mutationState)
{
    mMutationState = mutationState;
}


void TissueCell::SetSymmetricDivision()
{
    mSymmetricDivision = true;
}


bool TissueCell::DividesSymmetrically()
{ 
	return mSymmetricDivision;
}


void TissueCell::SetLogged()
{
    mIsLogged = true;
}


bool TissueCell::IsLogged()
{
    return mIsLogged;
}


/**
 * The TissueCell ready to divide method
 *
 */
bool TissueCell::ReadyToDivide()
{
    assert(!IsDead());
    if (mUndergoingApoptosis || mCellType==NECROTIC)
    {
        return false;
    }    
    
    mCanDivide = mpCellCycleModel->ReadyToDivide();

    return mCanDivide;
}


void TissueCell::StartApoptosis()
{
    assert(!IsDead());
    
    if (mUndergoingApoptosis)
    {
        EXCEPTION("StartApoptosis() called when already undergoing apoptosis");
    }
    mUndergoingApoptosis = true;
    
    CancerParameters *p_params = CancerParameters::Instance();
    
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    
    mDeathTime = p_simulation_time->GetDimensionalisedTime() + p_params->GetApoptosisTime();
}


bool TissueCell::HasApoptosisBegun() const
{
    return mUndergoingApoptosis;
}


double TissueCell::TimeUntilDeath() const
{
    if (!mUndergoingApoptosis)
    {
        EXCEPTION("Shouldn't be checking time until apoptosis as it isn't undergoing apoptosis");
    }
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    return mDeathTime - p_simulation_time->GetDimensionalisedTime();
}


bool TissueCell::IsDead() const
{
    SimulationTime *p_simulation_time = SimulationTime::Instance();

    return ( mIsDead || ( (mUndergoingApoptosis) && (p_simulation_time->GetDimensionalisedTime() >= mDeathTime)) );
}


void TissueCell::Kill()
{
    mIsDead = true;
}


void TissueCell::SetAncestor(unsigned ancestorIndex)
{
    mAncestor = ancestorIndex;
}


unsigned TissueCell::GetAncestor() const
{
    return mAncestor;   
}


TissueCell TissueCell::Divide()
{
    assert(!IsDead());
    
    // Copy this cell and give new one relevant attributes
    assert(mCanDivide);
    mCanDivide = false;
    
        
    if (mSymmetricDivision)
    {
        // Cell goes back to age zero, and cell type is possibly reset
        mpCellCycleModel->ResetModel();         
                
        TissueCell new_cell = TissueCell( 
            GetCellType(),  
            mMutationState, 
            mpCellCycleModel->CreateDaughterCellCycleModel()); 
                                         
        new_cell.GetCellCycleModel()->SetGeneration(1);
        new_cell.SetSymmetricDivision();
        new_cell.SetAncestor(GetAncestor());
        
        return new_cell;
    }
    else
    {
        mpCellCycleModel->SetGeneration(mpCellCycleModel->GetGeneration()+1);
        
        std::vector<CellType> new_cell_types = mpCellCycleModel->GetNewCellTypes();            
        mCellType = new_cell_types[0]; 
        
        // Cell goes back to age zero
        mpCellCycleModel->ResetModel();
        
        TissueCell new_cell = TissueCell( 
            new_cell_types[1],  
            mMutationState, 
            mpCellCycleModel->CreateDaughterCellCycleModel()); 
        
        assert(new_cell.GetCellCycleModel()->GetGeneration()==mpCellCycleModel->GetGeneration());
        mpCellCycleModel->SetMotherGeneration();
        new_cell.SetAncestor(GetAncestor());
        
        return new_cell;
    }        
}
