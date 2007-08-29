#ifndef _ODESECONDORDER_HPP_
#define _ODESECONDORDER_HPP_
#include "AbstractOdeSystem.hpp"


class OdeSecondOrder : public AbstractOdeSystem
{
public :
    OdeSecondOrder() : AbstractOdeSystem(2) // 2 here is the number of unknowns
    {
        mInitialConditions.push_back(0.0);
        mInitialConditions.push_back(1.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] =  rY[1];
        rDY[1] = -rY[0];
    }
};

#endif //_ODESECONDORDER_HPP_
