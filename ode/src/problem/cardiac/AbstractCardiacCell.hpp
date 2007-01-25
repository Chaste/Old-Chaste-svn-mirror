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
    const unsigned mVoltageIndex;  /**< The index of the voltage within our state variable vector */
    AbstractIvpOdeSolver *mpOdeSolver;   /**< Pointer to the solver used to simulate currents for this cell. */
    double mDt;
    AbstractStimulusFunction* mpIntracellularStimulus;
    AbstractStimulusFunction* mpExtracellularStimulus;
    
    // flag set to true if ComputeExceptVoltage is called
    bool mSetVoltageDerivativeToZero;
    
public:

    AbstractCardiacCell(AbstractIvpOdeSolver *pOdeSolver,
                        unsigned int numberOfStateVariables,
                        unsigned int voltageIndex,
                        double dt,
                        AbstractStimulusFunction* intracellularStimulus,
                        AbstractStimulusFunction* extracellularStimulus = NULL);
    
    virtual ~AbstractCardiacCell();
    
    /**
     * Initialise the cell:
     *   set our state variables to the initial conditions,
     *   set model parameters to their default values.
     */
    virtual void Init();
    
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep dt.
     */
    virtual OdeSolution Compute(double tStart, double tEnd);
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep dt, but does not update the voltage.
     */
    virtual OdeSolution ComputeExceptVoltage(double tStart, double tEnd);
    
    
    /**
     * Computes the total current flowing through the cell membrane, using the current
     * values of the state variables.
     * 
     * \todo does any cell model need to know the current time as well?
     */
    virtual double GetIIonic() = 0;
    
    void SetVoltage(double voltage);
    
    double GetVoltage();
    
    void SetStimulusFunction(AbstractStimulusFunction *stimulus);
    
    double GetStimulus(double time);
    
    void SetIntracellularStimulusFunction(AbstractStimulusFunction *stimulus);
    
    double GetIntracellularStimulus(double time);
    
    void SetExtracellularStimulusFunction(AbstractStimulusFunction *stimulus);
    
    double GetExtracellularStimulus(double time);
    
    bool HasExtracellularStimulus();
};

#endif /*ABSTRACTCARDIACCELL_HPP_*/
