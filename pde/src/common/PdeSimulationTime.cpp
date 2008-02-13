#include "PdeSimulationTime.hpp"

double PdeSimulationTime::mTime;

void PdeSimulationTime::SetTime(double time)
{
    mTime=time;
}

double PdeSimulationTime::GetTime()
{
    return mTime;
}
