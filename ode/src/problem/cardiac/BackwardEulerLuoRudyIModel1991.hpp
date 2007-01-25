#ifndef _BACKWARDEULERLUORUDYIMODEL1991_HPP_
#define _BACKWARDEULERLUORUDYIMODEL1991_HPP_

#include "AbstractStimulusFunction.hpp"
#include "AbstractBackwardEulerCardiacCell.hpp"

#include <vector>

/**
 * This class sets up the Luo-Rudy I 1991 system of equations, and solves them
 * using a decoupled backward Euler approach.
 */
class BackwardEulerLuoRudyIModel1991 : public AbstractBackwardEulerCardiacCell<1>
{
private:
    // Constants for the LuoRudyIModel1991OdeSystem model
    double background_current_E_b;
    double background_current_g_b;
    double fast_sodium_current_E_Na;
    double fast_sodium_current_g_Na;
    double ionic_concentrations_Ki;
    double ionic_concentrations_Ko;
    double ionic_concentrations_Nai;
    double ionic_concentrations_Nao;
    double membrane_C;
    double membrane_F;
    double membrane_R;
    double membrane_T;
    double plateau_potassium_current_g_Kp;
    double time_dependent_potassium_current_PR_NaK;
    
public:
    // Constructor
    BackwardEulerLuoRudyIModel1991(double dt,
                                   AbstractStimulusFunction *pIntracellularStimulus,
                                   AbstractStimulusFunction *pExtracellularStimulus = NULL);
                               
    // Destructor
    ~BackwardEulerLuoRudyIModel1991();
    
    void Init();
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep mDt.
     */
    OdeSolution Compute(double tStart, double tEnd);
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep mDt, but does not update the voltage.
     */
    OdeSolution ComputeExceptVoltage(double tStart, double tEnd);
    
protected:
    /**
     * Compute the values of all state variables except the voltage for one 
     * timestep.
     */
    void ComputeExceptVoltage(double tStart);

public:
    void ComputeResidual(const double rCurrentGuess[1], double rResidual[1]);
    void ComputeJacobian(const double rCurrentGuess[1], double rJacobian[1][1]);

    /**
     * Compute the ionic current at the current instant in time
     * (i.e. using the current values of the state variables).
     */
    double GetIIonic();
    
    /**
     *  Check that none of the gating variables have gone out of range. Throws an
     *  Exception if any have.
     */
    void VerifyGatingVariables();
};

#endif //
