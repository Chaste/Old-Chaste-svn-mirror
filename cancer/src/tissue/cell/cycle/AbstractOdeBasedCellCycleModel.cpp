#include "AbstractOdeBasedCellCycleModel.hpp"

AbstractOdeBasedCellCycleModel::AbstractOdeBasedCellCycleModel(double lastTime)
        : mpOdeSystem(NULL),
          mLastTime(lastTime),
          mDivideTime(lastTime),
          mReadyToDivide(false),
          mFinishedRunningOdes(false)
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


bool AbstractOdeBasedCellCycleModel::ReadyToDivide()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    
    double current_time = SimulationTime::Instance()->GetDimensionalisedTime();
    
//    we want to start solving the ODEs at the end of M phase - could possibly 
//    hijack mLastTime to do this
    
    if (current_time>mLastTime)
    {
        if (!mFinishedRunningOdes)
        {   
            mFinishedRunningOdes = SolveOdeToTime(current_time);
            
            for (unsigned i=0 ; i<mpOdeSystem->GetNumberOfStateVariables() ; i++)
            {
                if (mpOdeSystem->rGetStateVariables()[i]<0)
                {
                    #define COVERAGE_IGNORE
                    std::cout << "Protein["<< i <<"] = "<< mpOdeSystem->rGetStateVariables()[i] << "\n";
                    EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the CellCycleModel numerical method is probably unstable.");
                    #undef COVERAGE_IGNORE
                }
            }
            
//            mCurrentCellCyclePhase = G_ONE;
            
            if (mFinishedRunningOdes)
            {
//              mCurrentCellCyclePhase = S;  
                mDivideTime = GetOdeStopTime() + GetSG2Duration() + GetMDuration();
//              need to do some clever business here - instead of a divide time, we should
//              get back the separate pahse durations for S and G2
//              (need to do this in each concrete instance)
                if (current_time >= mDivideTime)
                {
                    mReadyToDivide = true;
                }
            }
            mLastTime = current_time;   // This is the last time the ODEs were evaluated.
        }
        else
        {   // ODE model finished, just increasing time until division...
            if (current_time >= mDivideTime)
            {
                mReadyToDivide = true;
            }
        }
    }
    return mReadyToDivide;
}


void AbstractOdeBasedCellCycleModel::ResetModel()
{
    assert(mReadyToDivide);
    mBirthTime = mDivideTime;
    mLastTime = mDivideTime;
    mFinishedRunningOdes = false;
    mReadyToDivide = false;
}
    
    
double AbstractOdeBasedCellCycleModel::GetSG2Duration()
{   
    // overridden in subclass StochasticWntCellCycleModel
    return CancerParameters::Instance()->GetSDuration()+CancerParameters::Instance()->GetG2Duration();
}   

double AbstractOdeBasedCellCycleModel::GetMDuration()
{   
    // overridden in subclass TysonNovak
    return CancerParameters::Instance()->GetMDuration();
}   
