#ifndef ABSTRACTCARDIACCELL_HPP_
#define ABSTRACTCARDIACCELL_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This is the base class for cardiac cell models. 
 */
class AbstractCardiacCell : public AbstractOdeSystem
{

protected:
    unsigned int mVoltageIndex;  /**< The index of the voltage within our state variable vector */ 
    AbstractIvpOdeSolver *mpOdeSolver;   /**< Pointer to the solver used to simulate currents for this cell. */
    double mDt;
    AbstractStimulusFunction* mpStimulus;
    
public:
    
    AbstractCardiacCell(AbstractIvpOdeSolver *pOdeSolver,
                        unsigned int numberOfStateVariables,
                        unsigned int voltageIndex, double dt)
        : AbstractOdeSystem(numberOfStateVariables)
    {
        mpOdeSolver = pOdeSolver;
        assert(voltageIndex < mNumberOfStateVariables);
        mVoltageIndex = voltageIndex;
        assert(dt>0);
        mDt=dt;
    }

    virtual ~AbstractCardiacCell()
    {
    }
    
    /**
     * Initialise the cell:
     *   set our state variables to the initial conditions,
     *   set model parameters to their default values.
     */
    virtual void Init()
    {
        rGetStateVariables() = mInitialConditions;
    }
    
    /**
     * Fudge to ensure that gating variables do not go out of bounds.
     */
    virtual void VerifyVariables() {}
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep dt.
     */
    virtual OdeSolution Compute(double tStart, double tEnd) 
    {
        return mpOdeSolver->Solve(this, rGetStateVariables(), tStart, tEnd, mDt, mDt);
    }
    
    /**
     * Computes the total current flowing through the cell membrane.
     */
    virtual double GetIIonic() = 0;
    
    void SetVoltage(double voltage)
    {
        rGetStateVariables()[mVoltageIndex] = voltage;
    }
    
    double GetVoltage()
    {
        return rGetStateVariables()[mVoltageIndex];
    }
    
    void SetStimulusFunction(AbstractStimulusFunction *stimulus)
    {
        mpStimulus = stimulus;
    }
     
    double GetStimulus(double time)
    {
        return mpStimulus->GetStimulus(time);
    }

    
    
    
    
    
    
    
    
    
    
};

#endif /*ABSTRACTCARDIACCELL_HPP_*/
