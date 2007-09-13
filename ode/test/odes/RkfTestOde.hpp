/*
 * Concrete class
 */
#ifndef _RKFTESTODE_HPP_
#define _RKFTESTODE_HPP_
#include "AbstractOdeSystem.hpp"

/*
 * Analytic solution to this, for y(0) = 0.5, is y = (t+1)^2 - 0.5*exp(t).
 */
class RkfTestOde : public AbstractOdeSystem
{
public :

    RkfTestOde() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mInitialConditions.push_back(0.5);
        mStateVariables = mInitialConditions;
    }
    
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0] - time*time + 1.0;
    }
};

#endif //_RKFTESTODE_HPP_
