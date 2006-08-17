#include "BackwardEulerIvpOdeSolverEvents.hpp"
#include <iostream>

bool BackwardEulerIvpOdeSolverEvents::CalculateStoppingEvent(AbstractOdeSystem* pAbstractOdeSystem, std::vector<double> currentYValue, double time)
{
    int num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();
    
    std::vector<double> dy(num_equations);
    dy = pAbstractOdeSystem->EvaluateYDerivatives(time, currentYValue);
    
    bool stopping_criterion;
    if (currentYValue[0] < 0.1 && dy[0] < 0.0)
    {
        stopping_criterion = true;
        std::cout << time << "\n";
    }
    else
    {
        stopping_criterion = false;
    }
    return stopping_criterion;
}
