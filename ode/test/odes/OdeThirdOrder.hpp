/**
 * Concrete OdeThirdOrder class
 */
#ifndef _ODETHIRDORDER_HPP
#define _ODETHIRDORDER_HPP
#include "AbstractOdeSystem.hpp"


class OdeThirdOrder : public AbstractOdeSystem
{
public :
    OdeThirdOrder()
            : AbstractOdeSystem(3) // 3 here is the number of unknowns
    {
        mInitialConditions.push_back(0.0);
        mInitialConditions.push_back(1.0);
        mInitialConditions.push_back(0.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]-rY[1]+rY[2];
        rDY[1]=rY[1]-rY[2];
        rDY[2]=2*rY[1]-rY[2];
    }
};

#endif //_ODETHIRDORDER_HPP
