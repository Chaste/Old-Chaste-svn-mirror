/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef TIMESTEPPER_HPP_
#define TIMESTEPPER_HPP_

/**
 * A helper class that provides a robust way to deal with time loops.
 *
 * An incremented integer counter is used to calculate the current time
 * and ensure that the correct end time.
 */
class TimeStepper
{
public:
    /**
     * Create a new time stepper over some simulation interval.
     * Time units are not specified, but all parameters should have consistent units.
     *
     * @param startTime  the start of the interval
     * @param endTime  the end of the interval
     * @param dt  the time step to use.
     * @param enforceConstantTimeStep If this is true the stepper estimates whether non-constant
     * timesteps will be used and quits if so.
     */
    TimeStepper(double startTime, double endTime, double dt, bool enforceConstantTimeStep=false);

    /**
     * Step forward one step in time and update member variables.
     */
    void AdvanceOneTimeStep();

    /**
     * Get the time
     */
    double GetTime() const;

    /**
     * Get the value time will take at the next step
     */
    double GetNextTime() const;

    /**
     * Get the size of the next time step which will be taken.
     * GetNextTimeStep() == GetNextTime() - GetTime()
     */
    double GetNextTimeStep();

    /**
     * True when GetTime ==  endTime
     */
    bool IsTimeAtEnd() const;

    /**
     * Estimate number of time steps, which may be an overestimate
     * Used to reserve memory for writing intermediate solutions.
     */
    unsigned EstimateTimeSteps() const;

    /**
     * The number of time AdvanceOneTimeStep() has been called.
     */
    unsigned GetTimeStepsElapsed() const;

private:
    friend class TestTimeStepper;

    /** The start time. */
    double mStart;

    /** The end time. */
    double mEnd;

    /** The size of time step. */
    double mDt;

    /** The number of time steps elapsed. */
    unsigned mTimeStep;

    /** The current time. */
    double mTime;

    /** What the time will be after the next time step. */
    double mNextTime;

    /** Compute what the time will be after the next time step. */
    double CalculateNextTime() const;
};

#endif /*TIMESTEPPER_HPP_*/
