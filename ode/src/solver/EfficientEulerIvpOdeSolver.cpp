#include "EfficientEulerIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

//#include <iostream>
#include <vector>
#include <cassert>

/**
 * Solve the given ODE system over the given domain, returning the solution
 * at intervals of duration timeSampling.
 *
 */
OdeSolution EfficientEulerIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem,
                          std::vector<double>& rInitialYValues,
                          double startTime,
                          double endTime,
                          double timeStep,
                          double timeSampling)
{
    // assert the size of the rYValues vector is correct
    const unsigned num_equations = rInitialYValues.size();
    assert(num_equations==pAbstractOdeSystem->GetNumberOfStateVariables());
    // Assert that we actually have a time interval > 0 .
    assert(endTime > startTime);
    // Assert that we  have a timestep > 0 .
    assert(timeStep > 0.0);
    // Assert that we  have a timesampling > 0 and >= timestep
    assert(timeSampling >= timeStep);

    // Determine the number of time steps that will be required to solve the
    // ODE system (note that the current algorithm accounts for any potential
    // floating point error)
    int number_of_time_samples = 0;
    double current_time = startTime;
    while (current_time < endTime)
    {
        number_of_time_samples++;
        if (startTime+number_of_time_samples*timeSampling >= endTime)
        {
            current_time = endTime;
        }
        else
        {
            current_time = startTime+number_of_time_samples*timeSampling;
        }
    }

    // Create solution object
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(number_of_time_samples);
    solutions.rGetTimes().push_back(startTime);
    std::vector<std::vector<double> >& all_y_values =
        solutions.rGetSolutions();
    all_y_values.push_back(rInitialYValues);
    
    // Allocate working memory
    std::vector<double> current_y_values(num_equations),
        next_y_values(num_equations);
    current_y_values.assign(rInitialYValues.begin(),
                            rInitialYValues.end());

    // Solve the ODE system
    int time_sample_number = 0;
    current_time = startTime;
    double end_time;
    bool curr_is_curr = false;

    while ( (current_time < endTime) )
    {
        time_sample_number++;
        end_time = startTime + time_sample_number*timeSampling;
        
        if (end_time >= endTime)
        {
            end_time = endTime;
        }
    
        { // Block used to be alternate Solve method call
            double real_time_step = timeStep;
            int time_step_number = 0;
            double start_time = current_time;
            while ( (current_time < end_time) )
            {
                time_step_number++;
                double to_time = start_time + time_step_number*timeStep;
                curr_is_curr = not curr_is_curr;
                
                // Determine what the value time step should really be like
                if (to_time >= end_time)
                {
                    real_time_step = end_time-current_time;
                    to_time = end_time;
                }
        
                // Function that calls the appropriate one-step solver
                CalculateNextYValue(pAbstractOdeSystem,
                                    real_time_step,
                                    current_time,
                                    curr_is_curr ? current_y_values : next_y_values,
                                    curr_is_curr ? next_y_values : current_y_values);
                                               
                // Determine the new current time
                current_time = to_time;
            }
        }
        
        //current_time = to_time; // Done in inner loop now
        
        // Push back new time into the time solution vector
        solutions.rGetTimes().push_back(current_time);
        // And new solution vector
        all_y_values.push_back(curr_is_curr ? next_y_values : current_y_values);
    }

    return solutions;
}


/**
 * Solve the given ODE system over the given domain.
 * Update the ODE's state with the final solution.
 *
 */
void EfficientEulerIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem,
                                       std::vector<double>& rInitialYValues,
                                       double startTime,
                                       double endTime,
                                       double timeStep)
{
    // assert the size of the rYValues vector is correct
    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();
    assert(rInitialYValues.size()==num_equations);
    // Assert that we actually have a time interval > 0 .
    assert(endTime > startTime);
    // Assert that we  have a timestep > 0 .
    assert(timeStep > 0.0);
    
    // Allocate working memory
    std::vector<double> current_y_values(num_equations),
        next_y_values(num_equations);
    current_y_values.assign(rInitialYValues.begin(),
                            rInitialYValues.end());

    // Solve the ODE system
    int time_step_number = 0;
    double current_time = startTime;
    double to_time;
    double real_time_step = timeStep;
    bool curr_is_curr = false;

    while ( (current_time < endTime) )
    {
        time_step_number++;
        curr_is_curr = not curr_is_curr;
        to_time = startTime + time_step_number*timeStep;
    
        // Determine what the value time step should really be like
        if (to_time >= endTime)
        {
            real_time_step = endTime-current_time;
            to_time = endTime;
        }
           
        // Function that calls the appropriate one-step solver
        CalculateNextYValue(pAbstractOdeSystem,
                            real_time_step,
                            current_time,
                            curr_is_curr ? current_y_values : next_y_values,
                            curr_is_curr ? next_y_values : current_y_values);
                                          
        // Determine the new current time
        current_time = to_time;
    }
    
    // Update ODE state
    pAbstractOdeSystem->SetStateVariables(curr_is_curr ? next_y_values : current_y_values);
}


/**
 * Calculate the solution for the next time step.
 *
 * Returns the solution in the argument rNextYValues, which must have
 * been pre-allocated to the correct size.
 *
 * @param pAbstractOdeSystem  ODE system to solve
 * @param timeStep  the current time step size
 * @param time  the current time
 * @param rCurrentYValues  the current solution
 * @param rNextYValues  vector to contain the next solution
 */
void EfficientEulerIvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             const std::vector<double>& rCurrentYValues,
                             std::vector<double>& rNextYValues)
{
    const unsigned num_equations = rNextYValues.size();
    assert(pAbstractOdeSystem->GetNumberOfStateVariables() == num_equations);
    assert(rCurrentYValues.size() == num_equations);
    
    // Evaluate derivatives, re-using solution memory
    pAbstractOdeSystem->EvaluateYDerivatives(time, rCurrentYValues, rNextYValues);

    // Update solution
    for (unsigned i=0; i<num_equations; i++)
    {
        rNextYValues[i] = rCurrentYValues[i] + timeStep*rNextYValues[i];
    }
}
