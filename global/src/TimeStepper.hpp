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
     * @param dt  the time step to use.  This must divide the simulation interval.
     */
    TimeStepper(double startTime, double endTime, double dt);

    void AdvanceOneTimeStep();
    
    double GetTime() const;
    double GetNextTime() const;    
    double GetNextTimeStep() const;
    
    bool IsTimeAtEnd() const;

private:
    double mStart;
    double mEnd;
    double mDt;
    unsigned mTimeStep;
    double mTime;
    double mNextTime;
    
    double CalculateNextTime();
};

#endif /*TIMESTEPPER_HPP_*/
