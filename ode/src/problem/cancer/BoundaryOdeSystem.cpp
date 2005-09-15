/**
 * Concrete BoundaryOde class
 */ 
#include "BoundaryOdeSystem.hpp"

BoundaryOdeSystem::BoundaryOdeSystem() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors

}

std::vector<double> BoundaryOdeSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(1);
	yDerivatives[0] = mCellVelocityAtBoundary;
	
	return yDerivatives;
}
	
