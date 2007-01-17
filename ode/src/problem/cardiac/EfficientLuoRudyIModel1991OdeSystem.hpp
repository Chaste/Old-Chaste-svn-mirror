#ifndef _EFFICIENTLUORUDYIMODEL1991ODESYSTEM_HPP_
#define _EFFICIENTLUORUDYIMODEL1991ODESYSTEM_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class sets up the EfficientLuoRudyIModel1991OdeSystem system of equations.
 */
class EfficientLuoRudyIModel1991OdeSystem : public AbstractCardiacCell
{
private:
    // Current and voltage components (objects) of the EfficientLuoRudyIModel1991OdeSystem model
    //   AbstractStimulusFunction *mpStimulus;
    
    // Constants for the EfficientLuoRudyIModel1991OdeSystem model
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
    EfficientLuoRudyIModel1991OdeSystem(AbstractIvpOdeSolver *pSolver,
                                        double dt,
                                        AbstractStimulusFunction *pIntracellularStimulus,
                                        AbstractStimulusFunction *pExtracellularStimulus = NULL);
                               
    // Destructor
    ~EfficientLuoRudyIModel1991OdeSystem();
    
    void Init();
    
    // This method will compute the RHS of the EfficientLuoRudyIModel1991OdeSystem model
    void EvaluateYDerivatives(double time, const std::vector<double> &rY,
                              std::vector<double>& rDY);
    
    double GetIIonic();
};

#endif //
