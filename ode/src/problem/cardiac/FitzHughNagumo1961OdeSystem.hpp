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
    /** 
     *  Constants for the FitzHugh-Nagumo model
     */
    static const double mAlpha = -0.08; 
    static const double mGamma = 3.00;
    static const double mEpsilon = 0.005;    

public:
    // Constructor
    FitzHughNagumo1961OdeSystem(AbstractIvpOdeSolver *pOdeSolver,
                                double dt,
                                AbstractStimulusFunction *pIntracelullarStimulus,
                                AbstractStimulusFunction *pExtracelullarStimulus=NULL);
                                
    // Destructor
    ~FitzHughNagumo1961OdeSystem();
    
    // Compute the RHS of the FitHugh-Nagumo system of ODEs
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY);
    double GetIIonic();
};

#endif //_FITZHUGHNAGUMO1961ODESYSTEM_HPP_
