/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef SIMULATIONTIME_HPP_
#define SIMULATIONTIME_HPP_

#include <boost/serialization/access.hpp>

#include "Exception.hpp"

/**
 * Simulation time object stores the simulation time.
 * It uses the singleton pattern to provide a globally consistent time.
 *
 * Note that the start time, end time and number of time steps must
 * be set before time can be incremented and returned.
 *
 * You should generally use the calls
 * IncrementTimeOneStep() and GetTime() when using this class.
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
    double GetTime() const;
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
    /**
     * A pointer to the singleton instance of this class.
     */
    static SimulationTime* mpInstance;
    /**
     * The duration of the simulation (cancer time units are in hours).
     */
    double mDurationOfSimulation;
    /**
     * The total number of steps for this simualation.
     */
    unsigned mTotalTimeStepsInSimulation;
    /**
     * The number of time steps which have been taken to date.
     */
    unsigned mTimeStepsElapsed;
    /**
     * A flag allowing us to determine whether the simulation time is ready to
     * be used.
     */
    bool mEndTimeAndNumberOfTimeStepsSet;
    /**
     * The current time (in hours)
     */
    double mCurrentTime;
    /**
     * The time at which the simulation should stop
     */
    double mEndTime;
    /**
     * Stores the time at which the simulation started
     */
    double mStartTime;
    /**
     * A flag allowing us to determine whether the start time of the simulation
     * has been set.
     */
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
        archive & mCurrentTime;
        archive & mEndTime;
        archive & mStartTimeSet;
        archive & mStartTime;
    }
};

#endif /*SIMULATIONTIME_HPP_*/
