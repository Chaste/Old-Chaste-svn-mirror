#ifndef _HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_
#define _HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_


#include "AbstractStimulusFunction.hpp"
#include "AbstractCardiacCell.hpp"
#include <vector>

/**
 * This class sets up the HodgkinHuxleySquidAxon1952OriginalOdeSystem system of equations.
 */
class HodgkinHuxleySquidAxon1952OriginalOdeSystem : public AbstractCardiacCell
{
private:
    /* Paramters */

    /*< mS/cm2 */
    static const double leakage_current_g_L = 0.3;    
    /*< uF/cm2 */
    static const double membrane_Cm = 1.0;    
    /*< mV */
    static const double membrane_E_R = -75.0;   
    /*< mS/cm2 */
    static const double potassium_channel_g_K = 36.0;  
    /*< mS/cm2 */
    static const double sodium_channel_g_Na = 120.0;   
    
public:
    // Constructor
    HodgkinHuxleySquidAxon1952OriginalOdeSystem(AbstractIvpOdeSolver *pOdeSolver,
                                                double dt,
                                                AbstractStimulusFunction *pIntracellularStimulus,
                                                AbstractStimulusFunction *pExtracellularStimulus=NULL);
    // Destructor
    ~HodgkinHuxleySquidAxon1952OriginalOdeSystem();
        
    // This method will compute the RHS of the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY);
    double GetIIonic();
};

#endif //_HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_
