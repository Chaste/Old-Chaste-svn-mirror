#include "AbstractOdeSystem.hpp"

// AbstractOdeSystem.cpp

AbstractOdeSystem::AbstractOdeSystem(const int& rNumberOfEquations,const double& rTInit, double * rYInit)
{
		mNumberOfEquations = rNumberOfEquations;
		mTInit = rTInit;
		mYInit = rYInit;
}

AbstractOdeSystem::~AbstractOdeSystem()
{
}
