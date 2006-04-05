#ifndef ABSTRACTCARDIACCELL_HPP_
#define ABSTRACTCARDIACCELL_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include <vector>

/**
 * This is the base class for cardiac cell models. 
 */
class AbstractCardiacCell : public AbstractOdeSystem
{

protected:
    unsigned int mVoltageIndex;  /**< The index of the voltage within our state variable vector */ 
    AbstractIvpOdeSolver *mpOdeSolver;   /**< Pointer to the solver used to simulate currents for this cell. */

public:
    
    AbstractCardiacCell(AbstractIvpOdeSolver *pOdeSolver,
                        unsigned int numberOfStateVariables,
                        unsigned int voltageIndex)
        : AbstractOdeSystem(numberOfStateVariables)
    {
        mpOdeSolver = pOdeSolver;
        assert(voltageIndex < mNumberOfStateVariables);
        mVoltageIndex = voltageIndex;
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
        mStateVariables = mInitialConditions;
    }
    
    /**
     * Fudge to ensure that gating variables do not go out of bounds.
     */
    virtual void VerifyVariables(std::vector<double>& odeVars) {}
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep dt.
     */
    virtual OdeSolution Compute(double tStart, double tEnd, double dt)
    {
        return mpOdeSolver->Solve(this, tStart, tEnd, dt, mStateVariables);
    }
    
    /**
     * Computes the total current flowing through the cell membrane.
     */
    virtual double GetIIonic() = 0;
    
    void SetVoltage(double voltage)
    {
        mStateVariables[mVoltageIndex] = voltage;
    }
    
    double GetVoltage()
    {
        return mStateVariables[mVoltageIndex];
    }
    
};

#endif /*ABSTRACTCARDIACCELL_HPP_*/
