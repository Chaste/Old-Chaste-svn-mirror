#ifndef _FORWARDEULERLUORUDYIMODEL1991_HPP_
#define _FORWARDEULERLUORUDYIMODEL1991_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>


/**
 * This class sets up the Luo-Rudy I 1991 system of equations.
 */
class ForwardEulerLuoRudyIModel1991 : public AbstractCardiacCell
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
    ForwardEulerLuoRudyIModel1991(double dt,
                                  AbstractStimulusFunction *pIntracellularStimulus,
                                  AbstractStimulusFunction *pExtracellularStimulus = NULL);
                               
    // Destructor
    ~ForwardEulerLuoRudyIModel1991();
    
    void Init();
    
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
    

    /** Compute the values of all state variables except the voltage for one 
     *  timestep
     */
    void ComputeExceptVoltage(double tStart);

    // This method will compute the RHS of the LuoRudyIModel1991OdeSystem model
    std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY);
    
    double GetIIonic();
    
    /**
     *  Check that none of the gating variables have gone out of range. Throws an
     *  Exception if any have.
     */
    void VerifyGatingVariables();
};

#endif //
