/**
 * Concrete RungeKuttaFehlbergIvpOdeSolver class.
 */
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <cmath>
#include <cfloat>
//#include <iostream>
#include <vector>

const double smidge=1e-10;

/*
 * PROTECTED FUNCTIONS
 */
 
 /**
  * This algorithm implements the `rkf45' routine using adaptive time stepping to 
  * solve as fast as possible to a given accuracy.
  * 
  * Note the change in the use of the solutions vector to give results when they are evaluated.
  * 
  * 
  */
void RungeKuttaFehlbergIvpOdeSolver::InternalSolve(OdeSolution& rSolution,
                                                AbstractOdeSystem* pOdeSystem,
                                                std::vector<double>& rYValues,
                                                std::vector<double>& rWorkingMemory,
                                                double startTime,
                                                double endTime,
                                                double maxTimeStep,
                                                double minTimeStep,
                                                double tolerance,
                                                bool outputSolution)
{
    double current_time = startTime;
    double time_step = maxTimeStep;
    bool got_to_end = false;
    bool accurate_enough = false;
    unsigned number_of_time_steps = 0;
    
    if (outputSolution)
    {   // Write out ICs
        rSolution.rGetTimes().push_back(current_time);
        rSolution.rGetSolutions().push_back(rYValues);
    }
    
    // should never get here if this bool has been set to true;
    assert(!mStoppingEventOccured);
    while ( !got_to_end )
    {
        //std::cout << "New timestep\n" << std::flush;
        while (!accurate_enough)
        {
            accurate_enough = true; // assume it is OK until we check and find otherwise

            // Function that calls the appropriate one-step solver
            std::vector<double> error = CalculateNextYValue(pOdeSystem,
                                time_step,
                                current_time,
                                rYValues,
                                rWorkingMemory);
              
            // Find the maximum error in this vector.
            double max_error = -DBL_MAX;
            for (unsigned i=0; i<error.size() ; i++)
            {
                if (error[i] > max_error)
                {
                    max_error = error[i];   
                }
            }
            
            if (max_error > tolerance)
            {// Reject the step-size and do it again.
                accurate_enough = false;
                //std::cout << "Approximation rejected\n" << std::flush;
            }
            else
            {
                // step forward the time since step has now been made
                current_time = current_time + time_step;    
                //std::cout << "Approximation accepted with time step = "<< time_step << "\n" << std::flush;
                //std::cout << "max_error = " << max_error << " tolerance = " << tolerance << "\n" << std::flush;
                if (outputSolution)
                {   // Write out ICs
                    //std::cout << "In solver Time = " << current_time << " y = " << rWorkingMemory[0] << "\n" << std::flush;
                    rSolution.rGetTimes().push_back(current_time);
                    rSolution.rGetSolutions().push_back(rWorkingMemory);
                    number_of_time_steps++;
                }
            }
            
            // Set a new step size based on the accuracy here
            AdjustStepSize(time_step, max_error, tolerance, maxTimeStep, minTimeStep);
        }
        
        // For the next timestep check the step doesn't go past the end...
        if (current_time + time_step > endTime)
        {   // Allow a smaller timestep for the final step.
            time_step = endTime - current_time;
        }
        
        if ( pOdeSystem->CalculateStoppingEvent(current_time,
                                                rWorkingMemory) == true )
        {
            mStoppingTime = current_time;
            mStoppingEventOccured = true;
        }
        
        if (mStoppingEventOccured || current_time>=endTime)
        {
            got_to_end = true;   
        }
        
        // Approximation accepted - put it in rYValues
        rYValues.assign(rWorkingMemory.begin(), rWorkingMemory.end());
        accurate_enough = false; // for next loop.
        //std::cout << "Finished Time Step\n" << std::flush;
    }
    rSolution.SetNumberOfTimeSteps(number_of_time_steps);
}


/**
 * Solves a system of ODEs using the Runge Kutta Fehlberg Adaptive timestep Initial Value Problem Ordinary Differential Equation Solver
 *
 * To be used in the form:
 *
 * RungeKuttaFehlbergIvpOdeSolver mySolver
 *
 * OdeSolution solution=mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, MaximumTimeStep, any_old_double);
 *
 * See documentation for AbstractIvpOdeSolver::Solve()
 * 
 * @return error the difference between the 4th and 5th order approximations.
 */
std::vector<double> RungeKuttaFehlbergIvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                                  double timeStep,
                                                  double time,
                                                  std::vector<double>& currentYValues,
                                                  std::vector<double>& nextYValues)
{
    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();

    std::vector<double> error(num_equations);
    std::vector<double> k1(num_equations);
    std::vector<double> k2(num_equations);
    std::vector<double> k3(num_equations);
    std::vector<double> k4(num_equations);
    std::vector<double> k5(num_equations);
    std::vector<double> k6(num_equations);
    std::vector<double>& dy = nextYValues; // re-use memory (not that it makes much difference here!)
    std::vector<double> yk2(num_equations);
    std::vector<double> yk3(num_equations);
    std::vector<double> yk4(num_equations);
    std::vector<double> yk5(num_equations);
    std::vector<double> yk6(num_equations);
    
    pAbstractOdeSystem->EvaluateYDerivatives(time, currentYValues, dy);
    
    for (unsigned i=0;i<num_equations; i++)
    {
        k1[i] = timeStep*dy[i];
        yk2[i] = currentYValues[i] + 0.25*k1[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time + 0.25*timeStep, yk2, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        k2[i] = timeStep*dy[i];
        yk3[i] = currentYValues[i] + 0.09375*k1[i] + 0.28125*k2[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time + 0.375*timeStep, yk3, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        k3[i] = timeStep*dy[i];
        yk4[i] = currentYValues[i] + m1932o2197*k1[i] - m7200o2197*k2[i] + m7296o2197*k3[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time+m12o13*timeStep, yk4, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        k4[i] = timeStep*dy[i];
        yk5[i] = currentYValues[i] + m439o216*k1[i] - 8*k2[i] 
                    + m3680o513*k3[i]- m845o4104*k4[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep, yk5, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        k5[i] = timeStep*dy[i];
        yk6[i] = currentYValues[i] - m8o27*k1[i] + 2*k2[i] - m3544o2565*k3[i]
                        + m1859o4104*k4[i] - 0.275*k5[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep, yk6, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        k6[i]=timeStep*dy[i];
        error[i] = (1/timeStep)*fabs(m1o360*k1[i] - m128o4275*k3[i]
                    - m2197o75240*k4[i] + 0.02*k5[i]+ m2o55*k6[i]);
        nextYValues[i] = currentYValues[i] + m25o216*k1[i] + m1408o2565*k3[i]
                        + m2197o4104*k4[i] - 0.2*k5[i];
    }
    return error;
}

void RungeKuttaFehlbergIvpOdeSolver::AdjustStepSize(double& rCurrentStepSize,
                                const double& rError, 
                                const double& rTolerance, 
                                const double& rMaxTimeStep, 
                                const double& rMinTimeStep)
{
    // Work out scaling factor delta for the step size
    double delta = pow(rTolerance/(2.0*rError), 0.25);
    
    // Maximum adjustment is *0.1 or *4
    if (delta <= 0.1)
    {
        rCurrentStepSize *= 0.1;
    }
    else if (delta >= 4.0)
    {
        rCurrentStepSize *= 4.0;  
    }
    else
    {
        rCurrentStepSize *= delta;
    }
    
    if (rCurrentStepSize > rMaxTimeStep)
    {
        rCurrentStepSize = rMaxTimeStep;
    }
    
    if (rCurrentStepSize < rMinTimeStep)
    {
        EXCEPTION("RKF45 Solver: Ode needs a smaller timestep than the set minimum\n");   
    }
    
}     


/*
 * PUBLIC FUNCTIONS
 */

OdeSolution RungeKuttaFehlbergIvpOdeSolver::Solve(AbstractOdeSystem* pOdeSystem,
                                               std::vector<double>& rYValues,
                                               double startTime,
                                               double endTime,
                                               double timeStep,
                                               double tolerance)
{
    assert(rYValues.size()==pOdeSystem->GetNumberOfStateVariables());
    assert(endTime > startTime);
    assert(timeStep > 0.0);
    
    mStoppingEventOccured = false;
    if ( pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("(Solve with sampling) Stopping event is true for initial condition");
    }
    // Allocate working memory
    std::vector<double> working_memory(rYValues.size());
    // And solve...
    OdeSolution solutions;
    //solutions.SetNumberOfTimeSteps((unsigned)(10.0*(startTime-endTime)/timeStep));
    bool return_solution = true;
    InternalSolve(solutions, pOdeSystem, rYValues, working_memory, startTime, endTime, timeStep, 1e-4, tolerance, return_solution);
    return solutions;
}

void RungeKuttaFehlbergIvpOdeSolver::Solve(AbstractOdeSystem* pOdeSystem,
                                        std::vector<double>& rYValues,
                                        double startTime,
                                        double endTime,
                                        double timeStep)
{
    assert(rYValues.size()==pOdeSystem->GetNumberOfStateVariables());
    assert(endTime > startTime);
    assert(timeStep > 0.0);
    
    mStoppingEventOccured = false;
    if ( pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("(Solve without sampling) Stopping event is true for initial condition");
    }
    // Allocate working memory
    std::vector<double> working_memory(rYValues.size());
    // And solve...
    OdeSolution not_required_solution;
    bool return_solution = false;
    InternalSolve(not_required_solution, pOdeSystem, rYValues, working_memory, startTime, endTime, timeStep, 1e-4, 1e-5, return_solution);
}


