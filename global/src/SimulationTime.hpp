#ifndef SIMULATIONTIME_HPP_
#define SIMULATIONTIME_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
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
    double mCurrentDimensionalisedTime;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        //archive & mpInstance;
        if(mpInstance==NULL)
        {
            EXCEPTION("Trying to save or load an instance of simulation time\n when the simulation time object does not exist");   
        }
        archive & mpInstance->mDurationOfSimulation;
        archive & mpInstance->mTotalTimeStepsInSimulation;
        archive & mpInstance->mTimeStepsElapsed;
        archive & mpInstance->mEndTimeAndNumberOfTimeStepsSet;
        archive & mpInstance->mCurrentDimensionalisedTime;
    }
};

#endif /*SIMULATIONTIME_HPP_*/
