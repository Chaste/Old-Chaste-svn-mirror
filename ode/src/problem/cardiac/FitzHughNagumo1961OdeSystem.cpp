#include "FitzHughNagumo1961OdeSystem.hpp"
#include <cmath>

/**
 * Constructor
 */
FitzHughNagumo1961OdeSystem::FitzHughNagumo1961OdeSystem(AbstractIvpOdeSolver *pOdeSolver,
                                                         double dt,
                                                         AbstractStimulusFunction *pIntracellularStimulus,
                                                         AbstractStimulusFunction *pExtracellularStimulus
                                                        )
        : AbstractCardiacCell(pOdeSolver,2,0,dt,pIntracellularStimulus,pExtracellularStimulus)
{
    /*
     * State variable
     */
    mVariableNames.push_back("V");
    mVariableUnits.push_back("mV");
    mInitialConditions.push_back(0.0);
    
    mVariableNames.push_back("w");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.0);
    Init();
}


/**
 * Destructor
 */
FitzHughNagumo1961OdeSystem::~FitzHughNagumo1961OdeSystem(void)
{
    // Do nothing
}


void FitzHughNagumo1961OdeSystem::Init()
{
    AbstractCardiacCell::Init();
    // Initialize model constants
    mAlpha = -0.08; // Typical values between 0.10 and 0.15
    mGamma = 3.00;
    mEpsilon = 0.005;
}

/**
 * Returns a vector representing the RHS of the FitzHugh-Nagumo system of Odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @return std::vector<double> RHS of FitzHugh-Nagumo system of equations
 */
std::vector<double> FitzHughNagumo1961OdeSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
    double membrane_V = rY[0]; // v
    double recovery_variable = rY[1]; // w
    
    double i_stim = GetStimulus(time);
    
    
    // dV/dt
    double membrane_V_prime = membrane_V*(membrane_V-mAlpha)*(1-membrane_V)-recovery_variable+i_stim;
    // do not update voltage if the mSetVoltageDerivativeToZero flag has been set
    if (mSetVoltageDerivativeToZero)
    {
        membrane_V_prime = 0;
    }
    
    
    // dw/dt
    double recovery_variable_prime = mEpsilon*(membrane_V-mGamma*recovery_variable);
    
    std::vector<double> RHS;
    
    RHS.push_back(membrane_V_prime);
    RHS.push_back(recovery_variable_prime);
    
    return RHS;
}

/**
 * Set the stimulus function used by this cell.
 *
 * @param stimulus  The stimulus function to use.
 */
//void FitzHughNagumo1961OdeSystem::SetStimulusFunction(AbstractStimulusFunction *stimulus)
//{
//    mpStimulus = stimulus;
//}

/**
 * Return the value of our stimulus function at the given time.
 */
//double FitzHughNagumo1961OdeSystem::GetStimulus(double time)
//{
//    return mpStimulus->GetStimulus(time);
//}


double FitzHughNagumo1961OdeSystem::GetIIonic()
{
    double membrane_V = mStateVariables[mVoltageIndex];
    double recovery_variable = mStateVariables[1];
    double fake_ionic_current = membrane_V*(membrane_V-mAlpha)*(1-membrane_V)-recovery_variable;
    return fake_ionic_current;
}
