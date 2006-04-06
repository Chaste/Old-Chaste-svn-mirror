#ifndef _FITZHUGHNAGUMO1961ODESYSTEM_HPP_
#define _FITZHUGHNAGUMO1961ODESYSTEM_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * Represents the FitzHugh-Nagumo system of ODEs.
 */
class FitzHughNagumo1961OdeSystem : public AbstractCardiacCell
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
      FitzHughNagumo1961OdeSystem(AbstractIvpOdeSolver *pOdeSolver, 
                                  AbstractStimulusFunction *stimulus, double dt);
                                  
      // Destructor
      ~FitzHughNagumo1961OdeSystem();

      void Init();
      
      void SetStimulusFunction(AbstractStimulusFunction *stimulus);
      double GetStimulus(double time);
      
      // Compute the RHS of the FitHugh-Nagumo system of ODEs
      std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY);
      void VerifyVariables(std::vector<double>& odeVars);
      double GetIIonic();
};

#endif //_FITZHUGHNAGUMO1961ODESYSTEM_HPP_
