#ifndef SIMULATIONTIME_HPP_
#define SIMULATIONTIME_HPP_

/**
 * Simulation time object stores the simulation time, uses the
 * singleton pattern
 */
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
    SimulationTime(const SimulationTime&);
    SimulationTime& operator= (const SimulationTime&);
private:
    static SimulationTime* mpInstance;
    double mDurationOfSimulation;
    int mTotalTimeStepsInSimulation;
    int mTimeStepsElapsed;
    bool mEndTimeAndNumberOfTimeStepsSet;
};

#endif /*SIMULATIONTIME_HPP_*/
