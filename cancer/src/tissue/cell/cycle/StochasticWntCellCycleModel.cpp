#include "StochasticWntCellCycleModel.hpp"

void StochasticWntCellCycleModel::SetG2Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double mean = CancerParameters::Instance()->GetG2Duration();
    double standard_deviation = 0.9;
    mG2Duration = p_gen->NormalRandomDeviate(mean,standard_deviation);
}

void StochasticWntCellCycleModel::InitialiseDaughterCell()
{
    SetG2Duration();
}

void StochasticWntCellCycleModel::Initialise()
{
    WntCellCycleModel::Initialise();
    SetG2Duration();
}

void StochasticWntCellCycleModel::ResetForDivision()
{
    AbstractWntOdeBasedCellCycleModel::ResetForDivision();
    SetG2Duration();
}

double StochasticWntCellCycleModel::GetG2Duration()
{
    return mG2Duration;
}

AbstractCellCycleModel* StochasticWntCellCycleModel::CreateDaughterCellCycleModel()
{
    assert(mpCell!=NULL);
    // calls a cheeky version of the constructor which makes the new cell cycle model
    // the same age as the old one - not a copy at this time.
    return new StochasticWntCellCycleModel(mpOdeSystem, mpCell->GetMutationState(), mBirthTime, mLastTime, mFinishedRunningOdes, mReadyToDivide,mDivideTime, mGeneration, mG2Duration);
}
