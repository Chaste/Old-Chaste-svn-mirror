#include "AbstractOdeBasedCellCycleModel.hpp"

AbstractOdeBasedCellCycleModel::AbstractOdeBasedCellCycleModel(double lastTime)
        : mpOdeSystem(NULL),
          mLastTime(lastTime),
          mDivideTime(lastTime),
          mReadyToDivide(false)
{
    AbstractCellCycleModel::SetBirthTime(lastTime);
}

AbstractOdeBasedCellCycleModel::~AbstractOdeBasedCellCycleModel()
{
    if (mpOdeSystem!=NULL)
    {
        delete mpOdeSystem;   
    }
}

void AbstractOdeBasedCellCycleModel::SetBirthTime(double birthTime)
{
    AbstractCellCycleModel::SetBirthTime(birthTime);
    mLastTime = birthTime;
    mDivideTime = birthTime;
}

std::vector<double> AbstractOdeBasedCellCycleModel::GetProteinConcentrations() const
{
    assert(mpOdeSystem!=NULL);
    return mpOdeSystem->rGetStateVariables();
}

void AbstractOdeBasedCellCycleModel::SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations)
{
    assert(mpOdeSystem!=NULL);
    assert(proteinConcentrations.size()==mpOdeSystem->rGetStateVariables().size());
    mLastTime = lastTime;
    mpOdeSystem->SetStateVariables(proteinConcentrations);
}

