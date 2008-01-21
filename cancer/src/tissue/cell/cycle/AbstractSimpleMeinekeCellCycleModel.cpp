#include "AbstractSimpleMeinekeCellCycleModel.hpp"


void AbstractSimpleMeinekeCellCycleModel::ResetForDivision()
{
    if (mGeneration+1u > CancerParameters::Instance()->GetMaxTransitGenerations())
    {
        mpCell->SetCellType(DIFFERENTIATED);
    }
    AbstractSimpleCellCycleModel::ResetForDivision();
    if (mGeneration == 1)
    {
        mGeneration = 0;
    }
}


void AbstractSimpleMeinekeCellCycleModel::InitialiseDaughterCell()
{
    if (mGeneration == 0)
    {
        mGeneration = 1;
    }
    // Daughter cell is always a TRANSIT or DIFFERENTIATED
    mpCell->SetCellType(TRANSIT);   
    if (mGeneration > CancerParameters::Instance()->GetMaxTransitGenerations())
    {
        mpCell->SetCellType(DIFFERENTIATED);
    }
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}
