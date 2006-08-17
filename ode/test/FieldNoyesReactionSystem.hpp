/**
 * Concrete FieldNoyesReactionSystem class
 */
#ifndef FIELDNOYESREACTIONSYSTEM_HPP_
#define FIELDNOYESREACTIONSYSTEM_HPP_
#include "AbstractOdeSystem.hpp"


class FieldNoyesReactionSystem : public AbstractOdeSystem
{
public :
    FieldNoyesReactionSystem() : AbstractOdeSystem(3)
    {
        mInitialConditions.push_back(1.0);
        mInitialConditions.push_back(1.0);
        mInitialConditions.push_back(1.0);
    }
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY)
    {
        std::vector<double> y_derivatives(3);
        const double epsilon = 0.05;
        const double p = 6.7;
        const double q = 1e-4;
        const double f = 0.5;
        
        y_derivatives[0] = (rY[1] - rY[0]*rY[1] + rY[0]*(1 - q*rY[0]))/epsilon;
        y_derivatives[1] = -rY[1] - rY[0]*rY[1] + 2*f*rY[2];
        y_derivatives[2] = (rY[0] - rY[2])/p;
        return y_derivatives;
    }
    
};

#endif /*FIELDNOYESREACTIONSYSTEM_HPP_*/
