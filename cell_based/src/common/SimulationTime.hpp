/*

Copyright (C) University of Oxford, 2005-2010

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

#include "ChasteSerialization.hpp"

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

    /**
     * Return a pointer to the simulation time object.
     * The first time this is called the simulation time object is created.
     */
    static SimulationTime* Instance();

    /**
     * Sets the end time and the number of time steps.
     * This must be called after SetStartTime() but before using any other methods.
     *
     * @param endTime time at which to end this run of the simulation
     * @param totalTimeStepsInSimulation  the number of time steps into which the above will be divided
     */
    void SetEndTimeAndNumberOfTimeSteps(double endTime, unsigned totalTimeStepsInSimulation);

    /**
     * Reset method for the end time and the number of time steps, to run the simulation
     * further after a first initial run.
     *
     * @param rEndTime the new end time for this simulation (the simulation will run from
     *      the current time to this new end time, NOT from 0 to this end time)
     * @param rNumberOfTimeStepsInThisRun the number of time steps into which the next run is split
     */
    void ResetEndTimeAndNumberOfTimeSteps(const double& rEndTime, const unsigned& rNumberOfTimeStepsInThisRun);

    /**
     * Get the simulation time step, set in earlier calls.
     *
     * Warning: Use of this method may result in round errors; generally use GetTime() instead.
     *
     * @return time step for this run of the simulation
     */
    double GetTimeStep() const;

    /**
     * Increment the simulation time by one time step.
     *
     * GetTime() will return an updated current time after this call.
     */
    void IncrementTimeOneStep();

    /**
     * Get the number of time steps that have elapsed.
     *
     * @return number of time steps which have been taken
     */
    unsigned GetTimeStepsElapsed() const;

    /**
     * Get the simulation time (in hours), should not have rounding errors.
     *
     * @return simulation time
     */
    double GetTime() const;

    /**
     * Destroy the current SimulationTime instance.  The next call to
     * Instance will create a new instance, on which
     * SetEndTimeAndNumberOfTimeSteps must be called again to reset time.
     *
     * This method *must* be called before program exit, to avoid a memory
     * leak.
     */
    static void Destroy();

    /**
     * Allows lower classes to check whether the simulation time class has been set up before using it
     *
     * @return whether the start time of the simulation has been set.
     */
    bool IsStartTimeSetUp() const;

    /**
     * @return whether the simulation has finished.
     */
    bool IsFinished() const;

    /**
     * @return the total number of time steps to be taken in this run.
     */
    unsigned GetTotalNumberOfTimeSteps() const;

    /**
     * Set the start time of the simulation
     *
     * @param startTime the time at which the simulation begins (usually 0.0 hours)
     */
    void SetStartTime(double startTime);

protected:

    /**
     * Default simulation time constructor
     *
     * Sets up time, you must set the start time,
     * end time and number of time steps before using the object.
     */
    SimulationTime();

    /**
     * Copy constructor.
     */
    SimulationTime(const SimulationTime&);

    /**
     * Overloaded assignment operator.
     */
    SimulationTime& operator= (const SimulationTime&);

private:

    /**
     * A pointer to the singleton instance of this class.
     */
    static SimulationTime* mpInstance;

    /**
     * The duration of the simulation (time is measured in units of hours).
     */
    double mDurationOfSimulation;

    /**
     * The total number of steps for this simulation.
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

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialization of a SimulationTime object must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
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
