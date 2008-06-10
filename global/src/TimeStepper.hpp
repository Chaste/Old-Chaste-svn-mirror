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


#ifndef TIMESTEPPER_HPP_
#define TIMESTEPPER_HPP_

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
     */
    TimeStepper(double startTime, double endTime, double dt);

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
    double GetNextTimeStep() const;

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
    double mStart;
    double mEnd;
    double mDt;
    unsigned mTimeStep;
    double mTime;
    double mNextTime;

    double CalculateNextTime() const;
};

#endif /*TIMESTEPPER_HPP_*/
