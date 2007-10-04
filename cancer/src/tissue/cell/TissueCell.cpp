#include "TissueCell.hpp"
#include "CellTypes.hpp"
#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"


TissueCell::TissueCell(CellType cellType,
                                   CellMutationState mutationState,
                                   unsigned generation,
                                   AbstractCellCycleModel *pCellCycleModel)
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
    // Stem cells are the only ones with generation = 0
    //assert( (generation == 0) == (cellType == STEM) ); Not for Wnt cells
    mGeneration = generation;
    mCellType = cellType;
    mMutationState = mutationState;
    mCanDivide = false;
    mUndergoingApoptosis = false;
    mIsDead = false;
    mDeathTime = DBL_MAX; // this has to be initialised for archiving...
    mNodeIndex = (unsigned)(-1); // initialise to a silly value for archiving (avoid memory check error)
    mIsLogged = false;
    mpCellCycleModel->SetCell(this);
}

void TissueCell::CommonCopy(const TissueCell &other_cell)
{
    // Copy 'easy' data members
    mGeneration = other_cell.mGeneration;
    mCellType = other_cell.mCellType;
    mMutationState = other_cell.mMutationState;
    mCanDivide = other_cell.mCanDivide;
    mUndergoingApoptosis = other_cell.mUndergoingApoptosis;
    mIsDead = other_cell.mIsDead;

    mDeathTime = other_cell.mDeathTime;
    mNodeIndex = other_cell.mNodeIndex;
    // Copy cell cycle model
    // First create a new object
    mpCellCycleModel = other_cell.mpCellCycleModel->CreateCellCycleModel();
    // Then copy its state
    *mpCellCycleModel = *(other_cell.mpCellCycleModel);

    mIsLogged = other_cell.mIsLogged;

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

void TissueCell::UpdateCellCycleModel()
{
    mpCellCycleModel->Update();
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

unsigned TissueCell::GetGeneration() const
{
    return mGeneration;
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
    if (mUndergoingApoptosis)
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


TissueCell TissueCell::Divide()
{
    assert(!IsDead());
    
    //Copy this cell and give new one relevant attributes...
    assert(mCanDivide);
    mCanDivide = false;
    
    CancerParameters *p_params = CancerParameters::Instance();
    //std::cout<< "Divide time" << mpSimulationTime->GetDimensionalisedTime() << "\n" ;
    if (mCellType != STEM)
    {
        if (mGeneration < p_params->GetMaxTransitGenerations())
        {
            mGeneration++;
            mpCellCycleModel->ResetModel();// Cell goes back to age zero
            return TissueCell(TRANSIT, mMutationState, mGeneration,
                                    mpCellCycleModel->CreateCellCycleModel());
        }
        else
        {
            mGeneration++;
            mCellType = DIFFERENTIATED;
            mpCellCycleModel->ResetModel();// Cell goes back to age zero
            return TissueCell(DIFFERENTIATED, mMutationState, mGeneration,
                                    mpCellCycleModel->CreateCellCycleModel());
        }
    }
    else
    {
        mpCellCycleModel->ResetModel();// Cell goes back to age zero
        return TissueCell(TRANSIT, mMutationState, 1,
                                mpCellCycleModel->CreateCellCycleModel());
    }
    
}

