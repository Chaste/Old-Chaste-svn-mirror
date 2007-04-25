#include "TimeStepper.hpp"
#include "Exception.hpp"
#include <cmath>
#include <cfloat>

TimeStepper::TimeStepper(double startTime, double endTime, double dt)
    : mStart(startTime),
      mTimeStep(0),
      mTime(startTime)
{
    if (startTime >= endTime)
    {
        EXCEPTION("The simulation duration must be strictly positive");
    }
    
    mTotalSteps = (unsigned) round ((endTime - startTime)/dt);
    mDt = (endTime - startTime) / mTotalSteps;
    
    if (fabs(mDt-dt) > DBL_EPSILON*10)
    {
        EXCEPTION("Can't create a time stepper unless dt divides simulation interval");
    }
}

unsigned TimeStepper::GetTotalSteps() const
{
    return mTotalSteps;
}

void TimeStepper::AdvanceOneTimeStep()
{
    mTimeStep++;
    mTime = mStart + mTimeStep * mDt;

}

double TimeStepper::GetTime() const
{
    return mTime;
}
