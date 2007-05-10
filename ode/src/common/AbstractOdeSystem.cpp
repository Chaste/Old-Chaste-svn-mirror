#include "AbstractOdeSystem.hpp"

AbstractOdeSystem::AbstractOdeSystem(unsigned numberOfStateVariables)
{
    mNumberOfStateVariables = numberOfStateVariables;
    mUseAnalytic = false;
}

AbstractOdeSystem::~AbstractOdeSystem()
{}

unsigned AbstractOdeSystem::GetStateVariableNumberByName(const std::string name)
{
    unsigned var_number=0;
    while (var_number != mNumberOfStateVariables &&
           mVariableNames[var_number] != name)
    {
        var_number++;
    }
    if (var_number == mNumberOfStateVariables)
    {
        EXCEPTION("State variable does not exist");
    }
    return var_number;
}
