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
        // Constants for the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
        double leakage_current_g_L;
        double membrane_Cm;
        double membrane_E_R;
        double potassium_channel_g_K;
        double sodium_channel_g_Na;

    public:
        // Constructor
        HodgkinHuxleySquidAxon1952OriginalOdeSystem(AbstractIvpOdeSolver *pOdeSolver, 
                                                    AbstractStimulusFunction *pStimulus, double dt);
        // Destructor
        ~HodgkinHuxleySquidAxon1952OriginalOdeSystem();
      
        void Init();
      
        // This method will compute the RHS of the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
        std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY);
        double GetIIonic();
};

#endif //_HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_
