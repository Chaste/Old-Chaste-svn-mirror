#ifndef _FITZHUGHNAGUMO1961ODESYSTEM_HPP_
#define _FITZHUGHNAGUMO1961ODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * Represents the FitzHugh-Nagumo system of ODEs.
 */
class FitzHughNagumo1961OdeSystem : public AbstractOdeSystem
{
   private:
      // Stimulus to be applied to cell
      AbstractStimulusFunction *mpStimulus;

      // Constants for the FitzHugh-Nagumo model
      double mAlpha;
      double mGamma;
      double mEpsilon;
      
   public:
      // Constructor
      FitzHughNagumo1961OdeSystem(AbstractStimulusFunction *stimulus);
      // Destructor
      ~FitzHughNagumo1961OdeSystem();

      void SetStimulusFunction(AbstractStimulusFunction *stimulus);
      double GetStimulus(double time);
      
      // Compute the RHS of the FitHugh-Nagumo system of ODEs
      std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY);
};

#endif //_FITZHUGHNAGUMO1961ODESYSTEM_HPP_
