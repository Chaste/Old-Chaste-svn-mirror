#include "AbstractOdeSystem.hpp"


AbstractOdeSystem::AbstractOdeSystem(unsigned numberOfStateVariables)
{
    mNumberOfStateVariables = numberOfStateVariables;
    mUseAnalytic = false;
}

AbstractOdeSystem::~AbstractOdeSystem()
{}
