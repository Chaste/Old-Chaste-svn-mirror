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
    void SetEndTimeAndNumberOfTimeSteps(double, unsigned);
    double GetTimeStep();
    void IncrementTimeOneStep();
    unsigned GetTimeStepsElapsed();
    double GetDimensionalisedTime();
    static void Destroy();
    bool IsSimulationTimeSetUp();
protected:
    SimulationTime();
    SimulationTime(const SimulationTime&);
    SimulationTime& operator= (const SimulationTime&);
private:
    static SimulationTime* mpInstance;
    double mDurationOfSimulation;
    unsigned mTotalTimeStepsInSimulation;
    unsigned mTimeStepsElapsed;
    bool mEndTimeAndNumberOfTimeStepsSet;
};

#endif /*SIMULATIONTIME_HPP_*/
