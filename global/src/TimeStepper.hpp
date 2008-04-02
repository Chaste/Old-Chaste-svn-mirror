/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
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
