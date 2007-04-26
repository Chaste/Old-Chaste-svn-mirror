#include "TimeStepper.hpp"
#include "Exception.hpp"
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <assert.h>

const double smidge=1e-10;

TimeStepper::TimeStepper(double startTime, double endTime, double dt)
    : mStart(startTime),
      mEnd(endTime),
      mDt(dt),
      mTimeStep(0),
      mTime(startTime)
{
    // std::cout << startTime << " " << endTime << " " << dt << std::endl;
    if (startTime >= endTime)
    {
        EXCEPTION("The simulation duration must be strictly positive");
    }
    mNextTime=CalculateNextTime();
}

double TimeStepper::CalculateNextTime() const
{
    double next_time=mStart + (mTimeStep+1) * mDt;
    if ((next_time) + smidge*(mDt) >= mEnd)
    {
        next_time = mEnd;
    }
    return next_time;
}

void TimeStepper::AdvanceOneTimeStep()
{   
    mTimeStep++;
    mTime = mNextTime;
    mNextTime=CalculateNextTime();
}

double TimeStepper::GetTime() const
{
    return mTime;
}

double TimeStepper::GetNextTime() const
{
    return mNextTime;
}

double TimeStepper::GetNextTimeStep() const
{
    double dt=mDt;
    if (mNextTime == mEnd)
    {
        dt = mEnd - mTime;
    }
    return dt;
}

bool TimeStepper::IsTimeAtEnd() const
{
    return mTime >= mEnd;
}

unsigned TimeStepper::EstimateTimeSteps() const
{
    return (unsigned) ceil((mEnd - mStart)/mDt);
}

unsigned TimeStepper::GetTimeStepsElapsed() const
{
    return mTimeStep;
}
