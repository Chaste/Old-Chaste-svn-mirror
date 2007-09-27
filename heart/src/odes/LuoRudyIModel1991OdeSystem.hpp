#ifndef _LUORUDYIMODEL1991ODESYSTEM_HPP_
#define _LUORUDYIMODEL1991ODESYSTEM_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class sets up the LuoRudyIModel1991OdeSystem system of equations.
 */
class LuoRudyIModel1991OdeSystem : public AbstractCardiacCell
{
private:

    /*< Constants for the LuoRudyIModel1991OdeSystem model */
    static const double membrane_C = 1.0;
    static const double membrane_F = 96484.6;
    static const double membrane_R = 8314;
    static const double membrane_T = 310.0;
    static const double background_current_E_b = -59.87;
    static const double background_current_g_b = 0.03921;
    static const double fast_sodium_current_g_Na = 23.0;
    static const double ionic_concentrations_Ki = 145.0;
    static const double ionic_concentrations_Ko = 5.4;
    static const double ionic_concentrations_Nai = 18.0;
    static const double ionic_concentrations_Nao = 140.0;
    static const double plateau_potassium_current_g_Kp = 0.0183;
    static const double time_dependent_potassium_current_PR_NaK = 0.01833;

    /*< Another parameter, which is a function of the above */
    double fast_sodium_current_E_Na;

    
public:
    // Constructor
    LuoRudyIModel1991OdeSystem(AbstractIvpOdeSolver *pSolver,
                               double dt,
                               AbstractStimulusFunction *pIntracellularStimulus,
                               AbstractStimulusFunction *pExtracellularStimulus = NULL);
                               
    // Destructor
    ~LuoRudyIModel1991OdeSystem();
        
    // This method will compute the RHS of the LuoRudyIModel1991OdeSystem model
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    double GetIIonic();
    
    double GetIntracellularCalciumConcentration();
};

#endif // _LUORUDYIMODEL1991ODESYSTEM_HPP_
