/**
 * Concrete AdamsBashforthIvpOdeSolver class. Sub-class of AbstractIvpOdeSolver.hpp
*/
#include "AdamsBashforthIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"
#include "Exception.hpp"

#include <iostream>
#include <vector>
#include <cassert>

/**
 * Solves a system of ODEs using the Adams-Bashforth Method Initial Value Problem Ordinary Differential Equation Solver
 *
 * @param pAbstractOdeSystem points to the concrete ODE system to be solved
 * @param rYValues a standard vector specifying the intial condition
 * of each solution variable in the system
 * @param startTime the time at which the initial conditions are specified
 * @param endTime the time to which the system should be solved and the solution
 * returned
 * @param timeStep the time interval to be used by the solver
 * @param timeSampling the time interval for generating the solution
 *
 * @return OdeSolution is an object containing an integer of the number of
 * equations, a std::vector of times and a std::vector of std::vectors where
 * each of those vectors contains the solution for one variable of the ODE
 * system at those times
 *
 * To be used in the form:
 *
 * AdamsBashforthIvpOdeSolver solver
 * OdeSolution solution=solver.Solve(pAbstractOdeSystem, rYValues, startTime, endTime, timeStep, timeSampling);
 *
 */

OdeSolution AdamsBashforthIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem,
                                              std::vector<double>& rYValues,
                                              double startTime,
                                              double endTime,
                                              double timeStep,
                                              double timeSampling)
{
    unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();
    
    // Assert that the size of Initial Conditions vector = number of equations.
    assert(rYValues.size()==num_equations);
    
    // Assert that we actually have a time interval > 0 .
    assert(endTime > startTime);
    
    // Assert that we  have a timestep > 0 .
    assert(timeStep > 0.0);
    
    // Assert that we  have a timesampling > 0 and >= timestep
    assert(timeSampling >= timeStep);
    
    mStoppingEventOccured = false;
    if ( pAbstractOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("Stopping event is true for initial condition");
    }
    
    // Determine the number of time steps that will be required to solve the
    // ODE system (note that the current algorithm accounts for any potential
    // floating point error)
    unsigned number_of_time_samples;
    double current_time;
    
    number_of_time_samples = 0;
    
    current_time = startTime;
    
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
    // number_of_time_samples has now been correctly set
    
    
    // setup output, write initial solution
    OdeSolution solutions;
    
    solutions.SetNumberOfTimeSteps(number_of_time_samples);
    solutions.rGetSolutions().push_back(rYValues);
    solutions.rGetTimes().push_back(startTime);
    
    // compute next writing (sampling) time
    unsigned time_sampling_counter = 1;
    double next_printing_time = startTime + time_sampling_counter*timeSampling;
    time_sampling_counter++;
    
    
    // Determine the number of time steps and make sure that we have at least
    // 4 of them
    unsigned number_of_timesteps = 0;
    
    current_time = startTime;
    while (current_time < endTime)
    {
        number_of_timesteps++;
        
        if (startTime+number_of_timesteps*timeStep >= endTime)
        {
            current_time = endTime;
        }
        else
        {
            current_time = startTime+number_of_timesteps*timeStep;
        }
    }
    // number_of_timesteps has now been correctly set
    
    if (number_of_timesteps <= 4)
    {
        EXCEPTION("A multi-step solver needs at least 4 time steps.");
    }
    
    std::vector<double> dy(num_equations);
    
    // set up a store for the derivatives of y at the 3 previous timesteps
    std::vector<std::vector<double> > derivative_store;
    derivative_store.resize(3);
    derivative_store[0].resize(num_equations);
    derivative_store[1].resize(num_equations);
    derivative_store[2].resize(num_equations);
    
    std::vector<double> dy_rk4(num_equations);
    std::vector<double> k1(num_equations);
    std::vector<double> k2(num_equations);
    std::vector<double> k3(num_equations);
    std::vector<double> k4(num_equations);
    
    std::vector<double> yk2(num_equations);
    std::vector<double> yk3(num_equations);
    std::vector<double> yk4(num_equations);
    
    
    // reset current time
    current_time = startTime;
    
    
    // Apply RungeKutta4's method first three timesteps, in order to
    // maintain fourth order accuracy of Adams-Bashforth method
    for (unsigned time_index=1; time_index<=3; time_index++)
    {
        pAbstractOdeSystem->EvaluateYDerivatives(current_time,rYValues,dy);
        
        for (unsigned i=0;i<num_equations; i++)
        {
            k1[i] = timeStep*dy[i];
            yk2[i] = rYValues[i] + 0.5*k1[i];
        }
        
        pAbstractOdeSystem->EvaluateYDerivatives(current_time+0.5*timeStep,yk2,dy_rk4);
        
        for (unsigned i=0;i<num_equations; i++)
        {
            k2[i] = timeStep*dy_rk4[i];
            yk3[i] = rYValues[i] + 0.5*k2[i];
        }
        
        pAbstractOdeSystem->EvaluateYDerivatives(current_time+0.5*timeStep,yk3,dy_rk4);
        
        for (unsigned i=0;i<num_equations; i++)
        {
            k3[i] = timeStep*dy_rk4[i];
            yk4[i] = rYValues[i] + k3[i];
        }

        pAbstractOdeSystem->EvaluateYDerivatives(current_time+timeStep,yk4,dy_rk4);
        
        for (unsigned i=0;i<num_equations; i++)
        {
            k4[i] = timeStep*dy_rk4[i];
            rYValues[i] = rYValues[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
        }
        
        // update current time
        current_time = startTime + time_index*timeStep;
        
        // Update OdeSolution if at next printing time (if current_time is either
        // greater than the next printing time OR current time is within 0.01% less
        // than the next printing time)
        if ( (current_time - next_printing_time)/next_printing_time > -0.0001)
        {
            solutions.rGetSolutions().push_back(rYValues);
            solutions.rGetTimes().push_back(current_time);
            next_printing_time = startTime + time_sampling_counter*timeSampling;
            time_sampling_counter++;
        }
        
        derivative_store[time_index-1]=dy;
    }
    
    
    double real_timestep = timeStep;
    
    unsigned timestep_number = 3;
    
    
    // Apply Adams-Bashforth method
    while ((current_time < endTime) && (!mStoppingEventOccured))
    {
        // Determine what the value time step should really be like
        double next_time = startTime+(timestep_number+1)*timeStep;
        if ( next_time >= endTime)
        {
            real_timestep = endTime-current_time;
        }
        
        pAbstractOdeSystem->EvaluateYDerivatives(current_time,rYValues,dy);
        
        for (unsigned i=0;i<num_equations; i++)
        {
            rYValues[i] = rYValues[i] + (real_timestep/24.0)*(55.0*dy[i] - 59.0*derivative_store[2][i] + 37.0*derivative_store[1][i] - 9.0*derivative_store[0][i]);
        }
        
        // Determine the new current time
        timestep_number++;
        if (real_timestep < timeStep)
        {
            current_time = endTime;
        }
        else
        {
            current_time = startTime+timestep_number*timeStep;
        }
        
        // check to see if a stopping event has occured
        if ( pAbstractOdeSystem->CalculateStoppingEvent(current_time, rYValues) == true )
        {
            mStoppingTime = current_time;
            mStoppingEventOccured = true;
            
            solutions.SetNumberOfTimeSteps(solutions.rGetTimes().size());
            solutions.rGetSolutions().push_back(rYValues);
            solutions.rGetTimes().push_back(current_time);
        }
        
        
        
        // Update OdeSolution if at next printing time (if current_time is either
        // greater than the next printing time OR current time is within 0.01% less
        // than the next printing time OR this is the final iteration)
        if ( ( (current_time - next_printing_time)/next_printing_time > -0.0001) ||
             (timestep_number == number_of_timesteps) )
        {
            solutions.rGetSolutions().push_back(rYValues);
            solutions.rGetTimes().push_back(current_time);
            //compute next printing time
            next_printing_time = startTime + time_sampling_counter*timeSampling;
            if (next_printing_time > endTime)
            {
                // this isn't strictly necessary as the above 'if' makes sure
                // output is written for the final iteration.
                next_printing_time = endTime;
            }
            time_sampling_counter++;
        }
        
        derivative_store[0] = derivative_store[1];
        derivative_store[1] = derivative_store[2];
        derivative_store[2] = dy;
        
    }
    
    return solutions;
}

/**
 * This method is required, since it occurs in the abstract base class.
 * However, it only makes much sense for the one step solvers.
 * Hence here it just behaves as the other solve method, and discards the OdeSolution.
 */
void AdamsBashforthIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem,
                                       std::vector<double>& rYValues,
                                       double startTime,
                                       double endTime,
                                       double timeStep)
{
    OdeSolution solution = Solve(pAbstractOdeSystem, rYValues, startTime, endTime, timeStep, 1.1*(endTime-startTime));
    rYValues = solution.rGetSolutions()[solution.GetNumberOfTimeSteps()];
}
