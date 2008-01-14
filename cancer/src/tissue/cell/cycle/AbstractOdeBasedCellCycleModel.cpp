#include "AbstractOdeBasedCellCycleModel.hpp"

AbstractOdeBasedCellCycleModel::AbstractOdeBasedCellCycleModel(double lastTime)
        : mpOdeSystem(NULL),
          mLastTime(lastTime),
          mDivideTime(lastTime),
          mFinishedRunningOdes(false),
          mG2PhaseStartTime(DBL_MAX)
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


void AbstractOdeBasedCellCycleModel::UpdateCellCyclePhase()
{     
    assert(mpOdeSystem!=NULL);
    
    double current_time = SimulationTime::Instance()->GetDimensionalisedTime();
    
    // Update the phase from M to G1 when necessary
    if (mCurrentCellCyclePhase == M_PHASE)
    {
        double m_duration = GetMDuration();
        if (GetAge() >= m_duration)
        {
            mCurrentCellCyclePhase = G_ONE_PHASE;
            mLastTime = m_duration + mBirthTime;
        }
        else
        {
            // Still dividing; don't run ODEs
            return;
        }
    }
    
    if (current_time > mLastTime)
    {
        if (!mFinishedRunningOdes)
        {   
            // Update whether a stopping event has occurred 
            mFinishedRunningOdes = SolveOdeToTime(current_time);
            
            // Check no concentrations have gone negative
            for (unsigned i=0 ; i<mpOdeSystem->GetNumberOfStateVariables() ; i++)
            {
                if (mpOdeSystem->rGetStateVariables()[i]< 0)
                {
                    #define COVERAGE_IGNORE
                    std::cout << "Protein["<< i <<"] = "<< mpOdeSystem->rGetStateVariables()[i] << "\n";
                    EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the CellCycleModel numerical method is probably unstable.");
                    #undef COVERAGE_IGNORE
                }
            }
                        
            if (mFinishedRunningOdes)
            {
                // Update durations of each phase
                mG1Duration = GetOdeStopTime() - mBirthTime - GetMDuration();
                mG2PhaseStartTime = GetOdeStopTime() + GetSDuration();
                mDivideTime = mG2PhaseStartTime + GetG2Duration();

                // Update phase
                if (current_time >= mG2PhaseStartTime)
                {
                    mCurrentCellCyclePhase = G_TWO_PHASE;
                }
                else
                {
                    mCurrentCellCyclePhase = S_PHASE;
                }
            }
            mLastTime = current_time;   // This is the last time the ODEs were evaluated.
        }
        else
        {
            // ODE model finished, just increasing time until division...
            if (current_time >= mG2PhaseStartTime)
            {
                mCurrentCellCyclePhase = G_TWO_PHASE;
            }
        }
    }
}


void AbstractOdeBasedCellCycleModel::ResetModel()
{
    assert(mFinishedRunningOdes);
    AbstractCellCycleModel::ResetModel();
    mBirthTime = mDivideTime;
    mLastTime = mDivideTime;
    mFinishedRunningOdes = false;
    mG1Duration = DBL_MAX;
    mDivideTime = DBL_MAX;
}
