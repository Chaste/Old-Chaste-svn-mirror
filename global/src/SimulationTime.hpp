#ifndef SIMULATIONTIME_HPP_
#define SIMULATIONTIME_HPP_

#include <boost/serialization/access.hpp>

#include "Exception.hpp"

/**
 * Simulation time object stores the simulation time, uses the
 * singleton pattern
 */
class SimulationTime
{
public:
    static SimulationTime* Instance();
    void SetEndTimeAndNumberOfTimeSteps(double, unsigned);
    void ResetEndTimeAndNumberOfTimeSteps(const double&, const unsigned&);
    double GetTimeStep();
    void IncrementTimeOneStep();
    unsigned GetTimeStepsElapsed();
    double GetDimensionalisedTime();
    static void Destroy();
    bool IsStartTimeSetUp();
    bool IsFinished();
    unsigned GetTotalNumberOfTimeSteps();
    void SetStartTime(double currentTime);
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
    double mCurrentDimensionalisedTime;
    double mEndTime;
    double mTimeAtEndOfLastRun;
    bool mStartTimeSet;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        //archive & mpInstance;
        if(mpInstance==NULL)
        {
            EXCEPTION("Trying to save or load an instance of simulation time when the simulation time object does not exist");   
        }
        archive & mpInstance->mDurationOfSimulation;
        archive & mpInstance->mTotalTimeStepsInSimulation;
        archive & mpInstance->mTimeStepsElapsed;
        archive & mpInstance->mEndTimeAndNumberOfTimeStepsSet;
        archive & mpInstance->mCurrentDimensionalisedTime;
        archive & mpInstance->mEndTime;
        archive & mpInstance->mStartTimeSet;
        archive & mpInstance->mTimeAtEndOfLastRun;
    }
};

#endif /*SIMULATIONTIME_HPP_*/
