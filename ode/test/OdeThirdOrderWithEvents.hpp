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
        
        
        std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY)
        {
            std::vector<double> y_derivatives(GetNumberOfStateVariables());
            y_derivatives[0]=rY[0]-rY[1]+rY[2];
            y_derivatives[1]=rY[1]-rY[2];
            y_derivatives[2]=2*rY[1]-rY[2];
            return y_derivatives;
        }
        
        bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
        {
            return (rY[0]<-0.5);
        }
    };
    
#endif //_ODETHIRDORDERWITHEVENTS_HPP
