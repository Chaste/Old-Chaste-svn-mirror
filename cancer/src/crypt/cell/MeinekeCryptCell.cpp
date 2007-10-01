#include "MeinekeCryptCell.hpp"
#include "MeinekeCryptCellTypes.hpp"
#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"


MeinekeCryptCell::MeinekeCryptCell(CryptCellType cellType,
                                   CryptCellMutationState mutationState,
                                   unsigned generation,
                                   AbstractCellCycleModel *pCellCycleModel)
        : mpCellCycleModel(pCellCycleModel)
{
    if (SimulationTime::Instance()->IsStartTimeSetUp()==false)
    {
        EXCEPTION("MeinekeCryptCell is setting up a cell cycle model but SimulationTime has not been set up");
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

void MeinekeCryptCell::CommonCopy(const MeinekeCryptCell &other_cell)
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

MeinekeCryptCell::MeinekeCryptCell(const MeinekeCryptCell &other_cell)
{
    CommonCopy(other_cell);
}

MeinekeCryptCell& MeinekeCryptCell::operator=(const MeinekeCryptCell &other_cell)
{
    // In case this is self-assignment, don't delete the cell cycle model...
    AbstractCellCycleModel* temp = mpCellCycleModel;
    CommonCopy(other_cell);
    // ...until after we've copied it.
    delete temp;
    return *this;
}

MeinekeCryptCell::~MeinekeCryptCell()
{
    delete mpCellCycleModel;
}


void MeinekeCryptCell::SetBirthTime(double birthTime)
{
    mpCellCycleModel->SetBirthTime(birthTime);
}

void MeinekeCryptCell::SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel)
{
    if (mpCellCycleModel != pCellCycleModel)
    {
        delete mpCellCycleModel;
    }
    mpCellCycleModel = pCellCycleModel;
    mpCellCycleModel->SetCell(this);
}

AbstractCellCycleModel* MeinekeCryptCell::GetCellCycleModel() const
{
    return mpCellCycleModel;
}

void MeinekeCryptCell::InitialiseCellCycleModel()
{
    mpCellCycleModel->Initialise();
}

/**
 * Set the node at which this cell is positioned.
 * 
 * @param index Index of the node
 */
void MeinekeCryptCell::SetNodeIndex(unsigned index)
{
    mNodeIndex = index;
}

unsigned MeinekeCryptCell::GetNodeIndex() const
{
    return mNodeIndex;
}

double MeinekeCryptCell::GetAge() const
{
    return mpCellCycleModel->GetAge();
}

double MeinekeCryptCell::GetBirthTime() const
{
    return mpCellCycleModel->GetBirthTime();
}

unsigned MeinekeCryptCell::GetGeneration() const
{
    return mGeneration;
}

CryptCellType MeinekeCryptCell::GetCellType() const
{
    return mCellType;
}

CryptCellMutationState MeinekeCryptCell::GetMutationState() const
{
    return mMutationState;
}

void MeinekeCryptCell::SetCellType(CryptCellType cellType)
{
    mCellType = cellType;
}

void MeinekeCryptCell::SetMutationState(CryptCellMutationState mutationState)
{
    mMutationState = mutationState;
}

void MeinekeCryptCell::SetLogged()
{
    mIsLogged = true;
}

bool MeinekeCryptCell::IsLogged()
{
    return mIsLogged;
}
/**
 * The MeinekeCryptCell ready to divide method
 *
 */
bool MeinekeCryptCell::ReadyToDivide()
{
    assert(!IsDead());
    if (mUndergoingApoptosis)
    {
        return false;
    }    
    
    mCanDivide = mpCellCycleModel->ReadyToDivide();

    return mCanDivide;
}


void MeinekeCryptCell::StartApoptosis()
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


bool MeinekeCryptCell::HasApoptosisBegun() const
{
    return mUndergoingApoptosis;
}

double MeinekeCryptCell::TimeUntilDeath() const
{
    if (!mUndergoingApoptosis)
    {
        EXCEPTION("Shouldn't be checking time until apoptosis as it isn't undergoing apoptosis");
    }
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    return mDeathTime - p_simulation_time->GetDimensionalisedTime();
}

bool MeinekeCryptCell::IsDead() const
{
    SimulationTime *p_simulation_time = SimulationTime::Instance();

    return ( mIsDead || ( (mUndergoingApoptosis) && (p_simulation_time->GetDimensionalisedTime() >= mDeathTime)) );
}

void MeinekeCryptCell::Kill()
{
    mIsDead = true;
}


MeinekeCryptCell MeinekeCryptCell::Divide()
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
            return MeinekeCryptCell(TRANSIT, mMutationState, mGeneration,
                                    mpCellCycleModel->CreateCellCycleModel());
        }
        else
        {
            mGeneration++;
            mCellType = DIFFERENTIATED;
            mpCellCycleModel->ResetModel();// Cell goes back to age zero
            return MeinekeCryptCell(DIFFERENTIATED, mMutationState, mGeneration,
                                    mpCellCycleModel->CreateCellCycleModel());
        }
    }
    else
    {
        mpCellCycleModel->ResetModel();// Cell goes back to age zero
        return MeinekeCryptCell(TRANSIT, mMutationState, 1,
                                mpCellCycleModel->CreateCellCycleModel());
    }
    
}

