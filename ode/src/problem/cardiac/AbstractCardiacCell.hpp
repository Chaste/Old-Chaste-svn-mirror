#ifndef ABSTRACTCARDIACCELL_HPP_
#define ABSTRACTCARDIACCELL_HPP_

#include "AbstractOdeSystem.hpp"
#include <vector>

/**
 * This is the base class for cardiac cell models. 
 */
class AbstractCardiacCell : public AbstractOdeSystem
{
private:
    unsigned int mVoltageIndex;  /**< The index of the voltage within our state variable vector */

protected:
    
    AbstractIvpOdeSolver *mpOdeSolver;   /**< Pointer to the solver used to simulate currents for this cell. */

public:
    
    AbstractCardiacCell(AbstractIvpOdeSolver *pOdeSolver,
                        unsigned int numberOfStateVariables,
                        unsigned int voltageIndex)
        : AbstractOdeSystem(numberOfStateVariables)
    {
        mpOdeSolver = pOdeSolver;
        asssert(voltageIndex < mNumberOfStateVariables);
        mVoltageIndex = voltageIndex;
    }

    virtual ~AbstractCardiacCell()
    {
    }
    
    /**
     * Fudge to ensure that gating variables do not go out of bounds.
     */
    virtual void VerifyVariables(std::vector<double>& odeVars) {}
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd].
     */
    virtual void Compute(double tStart, double tEnd)=0;
    
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
