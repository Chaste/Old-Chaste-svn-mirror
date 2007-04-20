#include "AbstractCardiacCell.hpp"

#include <cassert>

AbstractCardiacCell::AbstractCardiacCell(AbstractIvpOdeSolver *pOdeSolver,
                                         unsigned numberOfStateVariables,
                                         unsigned voltageIndex,
                                         double dt,
                                         AbstractStimulusFunction* intracellularStimulus,
                                         AbstractStimulusFunction* extracellularStimulus)
        : AbstractOdeSystem(numberOfStateVariables),
        mVoltageIndex(voltageIndex)
{
    mpOdeSolver = pOdeSolver;
    
    assert(voltageIndex < mNumberOfStateVariables);
    
    assert(dt>0);
    mDt=dt;
    
    mpIntracellularStimulus = intracellularStimulus;
    mpExtracellularStimulus = extracellularStimulus;
    
    mSetVoltageDerivativeToZero = false;
}

AbstractCardiacCell::~AbstractCardiacCell()
{}

/**
 * Initialise the cell:
 *   set our state variables to the initial conditions (done here),
 *   set model parameters to their default values (done by subclasses).
 */
void AbstractCardiacCell::Init()
{
    rGetStateVariables() = mInitialConditions;
}


/**
 * Simulates this cell's behaviour between the time interval [tStart, tEnd],
 * with timestep mDt.
 */
OdeSolution AbstractCardiacCell::Compute(double tStart, double tEnd)
{
    return mpOdeSolver->Solve(this, rGetStateVariables(), tStart, tEnd, mDt, mDt);
}


/**
 * Simulates this cell's behaviour between the time interval [tStart, tEnd],
 * with timestep mDt, but does not update the voltage.
 */
OdeSolution AbstractCardiacCell::ComputeExceptVoltage(double tStart, double tEnd)
{
    double saved_voltage = GetVoltage();
    
    mSetVoltageDerivativeToZero = true;
    OdeSolution ode_solution=mpOdeSolver->Solve(this, rGetStateVariables(), tStart, tEnd, mDt, mDt);
    mSetVoltageDerivativeToZero = false;
    
    SetVoltage(saved_voltage);
    return ode_solution;
}

void AbstractCardiacCell::SetVoltage(double voltage)
{
    rGetStateVariables()[mVoltageIndex] = voltage;
}

double AbstractCardiacCell::GetVoltage()
{
    return rGetStateVariables()[mVoltageIndex];
}

void AbstractCardiacCell::SetStimulusFunction(AbstractStimulusFunction *stimulus)
{
    SetIntracellularStimulusFunction(stimulus);
}

double AbstractCardiacCell::GetStimulus(double time)
{
    return GetIntracellularStimulus(time);
}

void AbstractCardiacCell::SetIntracellularStimulusFunction(AbstractStimulusFunction *stimulus)
{
    mpIntracellularStimulus = stimulus;
}

double AbstractCardiacCell::GetIntracellularStimulus(double time)
{
    return mpIntracellularStimulus->GetStimulus(time);
}

void AbstractCardiacCell::SetExtracellularStimulusFunction(AbstractStimulusFunction *stimulus)
{
    mpExtracellularStimulus = stimulus;
}

double AbstractCardiacCell::GetExtracellularStimulus(double time)
{
    /// \todo should this be an Exception?
    assert (HasExtracellularStimulus());
    
    return mpExtracellularStimulus->GetStimulus(time);
}

bool AbstractCardiacCell::HasExtracellularStimulus()
{
    return mpExtracellularStimulus != NULL;
}
