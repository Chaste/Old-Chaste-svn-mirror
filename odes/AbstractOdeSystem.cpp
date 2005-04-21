#include "AbstractOdeSystem.hpp"

// AbstractOdeSystem.cpp

AbstractOdeSystem::AbstractOdeSystem(const int& rNumberOfEquations)
{
		mNumberOfEquations = rNumberOfEquations;
		//mTInit = rTInit;
		//mYInit = rYInit;
}

AbstractOdeSystem::~AbstractOdeSystem()
{
}
