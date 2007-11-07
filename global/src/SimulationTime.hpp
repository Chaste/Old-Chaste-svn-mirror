#ifndef SIMULATIONTIME_HPP_
#define SIMULATIONTIME_HPP_

#include <boost/serialization/access.hpp>

#include "Exception.hpp"

/**
 * Simulation time object stores the simulation time.
 * It uses the singleton pattern to provide a globally consistent time.
 */
class SimulationTime
{
public:
    static SimulationTime* Instance();
    void SetEndTimeAndNumberOfTimeSteps(double, unsigned);
    void ResetEndTimeAndNumberOfTimeSteps(const double&, const unsigned&);
    double GetTimeStep() const;
    void IncrementTimeOneStep();
    unsigned GetTimeStepsElapsed() const;
    double GetDimensionalisedTime() const;
    static void Destroy();
    bool IsStartTimeSetUp() const;
    bool IsFinished() const;
    unsigned GetTotalNumberOfTimeSteps() const;
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
    /**
     * Serialization of a SimulationTime object must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & mDurationOfSimulation;
        archive & mTotalTimeStepsInSimulation;
        archive & mTimeStepsElapsed;
        archive & mEndTimeAndNumberOfTimeStepsSet;
        archive & mCurrentDimensionalisedTime;
        archive & mEndTime;
        archive & mStartTimeSet;
        archive & mTimeAtEndOfLastRun;
    }
};

#endif /*SIMULATIONTIME_HPP_*/
