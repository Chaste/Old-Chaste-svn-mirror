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

std::string AbstractOdeSystem::DumpState(const std::string& message)
{
    std::stringstream res;
    res << message << "\nState:\n";
    const std::vector<double>& rY = rGetStateVariables();
    for (unsigned i=0; i<mVariableNames.size(); i++)
    {
        res << "\t" << mVariableNames[i] << ":" << rY[i] << "\n";
    }
    return res.str();
}
