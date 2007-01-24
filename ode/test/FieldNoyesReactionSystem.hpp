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
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        const double epsilon = 0.05;
        const double p = 6.7;
        const double q = 1e-4;
        const double f = 0.5;
        
        rDY[0] = (rY[1] - rY[0]*rY[1] + rY[0]*(1 - q*rY[0]))/epsilon;
        rDY[1] = -rY[1] - rY[0]*rY[1] + 2*f*rY[2];
        rDY[2] = (rY[0] - rY[2])/p;
    }
    
};

#endif /*FIELDNOYESREACTIONSYSTEM_HPP_*/
