#include "AbstractCellCycleModel.hpp"

AbstractCellCycleModel::~AbstractCellCycleModel()
{
    // Don't delete the cell - the cell deletes the cell cycle model 
    // when it is destroyed (not the following way round):
    // delete mpCell;
}

void AbstractCellCycleModel::SetCell(TissueCell* pCell)
{
    mpCell = pCell;
}

TissueCell* AbstractCellCycleModel::GetCell()
{
    assert(mpCell!=NULL);
    return mpCell;
}

void AbstractCellCycleModel::SetBirthTime(double birthTime)
{
    mBirthTime = birthTime;
}

double AbstractCellCycleModel::GetBirthTime() const
{
    return mBirthTime;
}

double AbstractCellCycleModel::GetAge()
{
    return SimulationTime::Instance()->GetDimensionalisedTime() - mBirthTime;
}

CellCyclePhase AbstractCellCycleModel::GetCurrentCellCyclePhase()
{
    return mCurrentCellCyclePhase;
}

#define COVERAGE_IGNORE
double AbstractCellCycleModel::GetMembraneBoundBetaCateninLevel()
{
    NEVER_REACHED;
    return 0.0;
}

double AbstractCellCycleModel::GetCytoplasmicBetaCateninLevel()
{
    NEVER_REACHED;
    return 0.0;
}

double AbstractCellCycleModel::GetNuclearBetaCateninLevel()
{
    NEVER_REACHED;
    return 0.0;
}
#undef COVERAGE_IGNORE

bool AbstractCellCycleModel::UsesBetaCat()
{
    return false;
}

void AbstractCellCycleModel::SetGeneration(unsigned generation)
{
    mGeneration = generation;
}

unsigned AbstractCellCycleModel::GetGeneration() const
{
    return mGeneration;
}

void AbstractCellCycleModel::ResetModel()
{
    assert(mReadyToDivide);
    mGeneration++;
    mCurrentCellCyclePhase = M_PHASE;
    mReadyToDivide = false;
}
    
double AbstractCellCycleModel::GetSDuration()
{
    return CancerParameters::Instance()->GetSDuration();
}   

double AbstractCellCycleModel::GetG1Duration()
{
    return mG1Duration;
}

double AbstractCellCycleModel::GetG2Duration()
{
    return CancerParameters::Instance()->GetG2Duration();
}   

double AbstractCellCycleModel::GetMDuration()
{
    return CancerParameters::Instance()->GetMDuration();
}   

bool AbstractCellCycleModel::ReadyToDivide()
{
    assert(mpCell != NULL);
    
    if (!mReadyToDivide)
    {
        UpdateCellCyclePhase();
        if ( (mCurrentCellCyclePhase != G_ZERO_PHASE) &&
             (GetAge() >= GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration()) )
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
}
