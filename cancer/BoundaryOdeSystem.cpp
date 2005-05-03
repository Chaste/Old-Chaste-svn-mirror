/**
 * Concrete BoundaryOde class
 */ 
#include "BoundaryOdeSystem.hpp"

BoundaryOdeSystem::BoundaryOdeSystem() : AbstractOdeSystem(1)
{
// Use AbstractOdeSystem constructors

}

std::vector<double> BoundaryOdeSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(mNumberOfEquations);
	yDerivatives[0] = mVelocityCellAtBoundary;
	
	return yDerivatives;
}
	
