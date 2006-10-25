#ifndef SIMULATIONTIME_HPP_
#define SIMULATIONTIME_HPP_

class SimulationTime
{
public:
    static SimulationTime* Instance();
    void SetEndTimeAndNumberOfTimeSteps(double, int);
    double GetTimeStep();
    void IncrementTimeOneStep();
    int GetTimeStepsElapsed();
    double GetDimensionalisedTime();
    static void Destroy();
protected:
    SimulationTime();
private:
    static SimulationTime* mpInstance;
    double mDurationOfSimulation;
    int mTotalTimeStepsInSimulation;
    int mTimeStepsElapsed;
    bool mEndTimeAndNumberOfTimeStepsSet;
};

#endif /*SIMULATIONTIME_HPP_*/
