#ifndef _LUORUDYIMODEL1991ODESYSTEM_HPP_
#define _LUORUDYIMODEL1991ODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class sets up the LuoRudyIModel1991OdeSystem system of equations.
 */
class LuoRudyIModel1991OdeSystem : public AbstractOdeSystem
{
   private:
      // Current and voltage components (objects) of the LuoRudyIModel1991OdeSystem model
      AbstractStimulusFunction *mpStimulus;

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
      LuoRudyIModel1991OdeSystem(AbstractStimulusFunction *stimulus);
      // Destructor
      ~LuoRudyIModel1991OdeSystem();

      // This method will compute the RHS of the LuoRudyIModel1991OdeSystem model
      std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY);
};

#endif //
