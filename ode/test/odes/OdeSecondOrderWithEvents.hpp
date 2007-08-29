#ifndef _ODESECONDORDERWITHEVENTS_HPP_
#define _ODESECONDORDERWITHEVENTS_HPP_
#include "AbstractOdeSystem.hpp"


class OdeSecondOrderWithEvents : public AbstractOdeSystem
{
public :
    OdeSecondOrderWithEvents()
            : AbstractOdeSystem(2)  // 2 here is the number of variables
    {
        // set initial conditions
        mInitialConditions.push_back(1.0);
        mInitialConditions.push_back(0.0);
    }
    
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] =  rY[1];
        rDY[1] = -rY[0];
    }
    
    
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {
        return (rY[0]<0);
    }
};

#endif //_ODESECONDORDERWITHEVENTS_HPP_
