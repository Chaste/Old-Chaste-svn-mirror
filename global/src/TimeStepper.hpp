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
    
    /**
     * Return the number of time steps that will be used in total.
     */
    unsigned GetTotalSteps() const;
    
    void AdvanceOneTimeStep();
    
    double GetTime() const;

private:
    double   mStart;
    unsigned mTotalSteps;
    unsigned mTimeStep;
    double mDt;
    double mTime;
};

#endif /*TIMESTEPPER_HPP_*/
