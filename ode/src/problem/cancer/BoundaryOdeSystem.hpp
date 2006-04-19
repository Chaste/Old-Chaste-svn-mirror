/**
 * Concrete BoundaryOde class
 */ 
#ifndef _BOUNDARYODESYSTEM_HPP_
#define _BOUNDARYODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"


class BoundaryOdeSystem : public AbstractOdeSystem
{
private:
    double mCellVelocityAtBoundary;

public:
    
    BoundaryOdeSystem();
    
    std::vector<double> BoundaryOdeSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY);


};

#endif //_BOUNDARYODESYSTEM_HPP_
