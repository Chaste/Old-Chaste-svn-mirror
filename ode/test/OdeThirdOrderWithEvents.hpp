#/**
  * Concrete OdeThirdOrder class with events
  */
#ifndef _ODETHIRDORDERWITHEVENTS_HPP
#define _ODETHIRDORDERWITHEVENTS_HPP
#include "AbstractOdeSystem.hpp"
            
            
class OdeThirdOrderWithEvents : public AbstractOdeSystem
{
public :
    OdeThirdOrderWithEvents()
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
    
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {
        return (rY[0]<-0.5);
    }
};
    
#endif //_ODETHIRDORDERWITHEVENTS_HPP
