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

std::vector<CellType> AbstractCellCycleModel::GetNewCellTypes(CellType cellType)
{
    std::vector<CellType> cell_types(2);
    if (cellType == STEM)
    {
        cell_types[0] = cellType;
        cell_types[1] = TRANSIT;
    }
    else
    {
        cell_types[0] = cellType;
        cell_types[1] = cellType;
    }
    return cell_types;
}

void AbstractCellCycleModel::SetMotherGeneration(CellType cellType)
{
    if (cellType == STEM)
    {
        mGeneration--;
    }
}

   
double AbstractCellCycleModel::GetSDuration()
{
    return CancerParameters::Instance()->GetSDuration();
}   
    
double AbstractCellCycleModel::GetG2Duration()
{
    return CancerParameters::Instance()->GetG2Duration();
}   

double AbstractCellCycleModel::GetMDuration()
{
    return CancerParameters::Instance()->GetMDuration();
}   

