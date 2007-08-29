#ifndef VANDERPOLODE_HPP_
#define VANDERPOLODE_HPP_

#include "AbstractOdeSystem.hpp"


class VanDerPolOde : public AbstractOdeSystem
{
public :
    VanDerPolOde() : AbstractOdeSystem(2)  // 2 here is the number of unknowns
    {
        mInitialConditions.push_back(10.0);
        mInitialConditions.push_back(10.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        double mu = 1;
        rDY[0]= rY[1] + mu*(rY[0] - rY[0]*rY[0]*rY[0]);
        rDY[1] = -rY[0];
    }
    
};


#endif /*VANDERPOLODE_HPP_*/
