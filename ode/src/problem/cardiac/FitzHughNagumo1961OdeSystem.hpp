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
    // Constants for the FitzHugh-Nagumo model
    double mAlpha;
    double mGamma;
    double mEpsilon;
    
public:
    // Constructor
    FitzHughNagumo1961OdeSystem(AbstractIvpOdeSolver *pOdeSolver,
                                double dt,
                                AbstractStimulusFunction *pIntracelullarStimulus,
                                AbstractStimulusFunction *pExtracelullarStimulus=NULL);
                                
    // Destructor
    ~FitzHughNagumo1961OdeSystem();
    
    void Init();
    
    // Compute the RHS of the FitHugh-Nagumo system of ODEs
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY);
    double GetIIonic();
};

#endif //_FITZHUGHNAGUMO1961ODESYSTEM_HPP_
