#ifndef ABSTRACTBACKWARDEULERCARDIACCELL_HPP_
#define ABSTRACTBACKWARDEULERCARDIACCELL_HPP_

#include "AbstractCardiacCell.hpp"
//#include "CardiacNewtonSolver.hpp"

#include <cassert>
#include <cmath>

/**
 * This is the base class for cardiac cells solved using a (decoupled) backward
 * Euler approach.
 *
 * The basic approach to solving such models is:
 *  * Update the transmembrane potential, either from solving an external PDE,
 *    or using a forward Euler step.
 *  * Update any gating variables (or similar) using a backward euler step.
 *    Suitable ODEs can be written in the form  du/dt = g(V) + h(V)*u.  The update
 *    expression is then  u_n = ( u_{n-1} + g(V_n)*dt ) / ( 1 - h(V_n)*dt ).
 *  * Update the remaining state variables using Newton's method to solve the
 *    nonlinear system  U_n - U_{n-1} = dt*F(U_n, V_n).
 *    The template parameter to the class specifies the size of this nonlinear system.
 */
template<unsigned SIZE>
class AbstractBackwardEulerCardiacCell : public AbstractCardiacCell
{
public:

    /**
     * Standard constructor for a cell.
     * 
     * @param numberOfStateVariables  the size of the ODE system
     * @param voltageIndex  the index of the variable representing the transmembrane
     *     potential within the state variable vector
     * @param dt  the timestep to use in solving the ODEs
     * @param intracellularStimulus  the intracellular stimulus function
     * @param extracellularStimulus  the extracellular stimulus function
     * 
     * Some notes for future reference:
     *  * We may want to remove the timestep from this class, and instead pass it to
     *    the Compute* methods, especially if variable timestepping is to be used.
     *  * It's a pity that inheriting from AbstractCardiacCell forces us to store a
     *    null pointer (for the unused ODE solver) in every instance.  We may want
     *    to revisit this design decision at a later date.
     */
    AbstractBackwardEulerCardiacCell(
        unsigned numberOfStateVariables,
        unsigned voltageIndex,
        double dt,
        AbstractStimulusFunction* intracellularStimulus,
        AbstractStimulusFunction* extracellularStimulus = NULL);

    virtual ~AbstractBackwardEulerCardiacCell();
    
    /**
     * Compute the residual of the nonlinear system portion of the cell model.
     * 
     * @param rCurrentGuess  the current guess for U_n
     * @param rResidual  to be filled in with the residual vector
     */
    virtual void ComputeResidual(const double rCurrentGuess[SIZE], double rResidual[SIZE])=0;
    
    /**
     * Compute the Jacobian matrix for the nonlinear system portion of the cell model.
     * 
     * @param rCurrentGuess  the current guess for U_n
     * @param rJacobian  to be filled in with the Jacobian matrix
     */
    virtual void ComputeJacobian(const double rCurrentGuess[SIZE], double rJacobian[SIZE][SIZE])=0;
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep mDt.  Uses a forward Euler step to update the transmembrane
     * potential at each timestep.
     * 
     * The length of the time interval must be a multiple of the timestep.
     * 
     * @return  the values of each state variable, at mDt intervals.
     */
    virtual OdeSolution Compute(double tStart, double tEnd);
    
    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep mDt.  The transmembrane potential is kept fixed throughout.
     * 
     * The length of the time interval must be a multiple of the timestep.
     */
    virtual void ComputeExceptVoltage(double tStart, double tEnd);
    
    /**
     *  Check that none of the gating variables have gone out of range. Throws an
     *  Exception if any have.
     */
    virtual void VerifyGatingVariables()=0;
    
private:
#define COVERAGE_IGNORE
    /**
     * This function should never be called - the cell class incorporates its own solver.
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
    {
        NEVER_REACHED;
    }
#undef COVERAGE_IGNORE
    
protected:
    /**
     * Compute the values of all state variables except the voltage, for one 
     * timestep from tStart.
     * 
     * This method must be provided by subclasses.
     */
    virtual void ComputeOneStepExceptVoltage(double tStart)=0;
    
    /**
     * Perform a forward Euler step to update the transmembrane potential.
     * 
     * This method must be provided by subclasses.
     */
    virtual void UpdateTransmembranePotential(double time)=0;
};

template <unsigned SIZE>
AbstractBackwardEulerCardiacCell<SIZE>::AbstractBackwardEulerCardiacCell(
    unsigned numberOfStateVariables,
    unsigned voltageIndex,
    double dt,
    AbstractStimulusFunction* intracellularStimulus,
    AbstractStimulusFunction* extracellularStimulus)
        : AbstractCardiacCell(NULL,
                              numberOfStateVariables,
                              voltageIndex,
                              dt,
                              intracellularStimulus,
                              extracellularStimulus)
{}

template <unsigned SIZE>
AbstractBackwardEulerCardiacCell<SIZE>::~AbstractBackwardEulerCardiacCell()
{}

template <unsigned SIZE>
OdeSolution AbstractBackwardEulerCardiacCell<SIZE>::Compute(double tStart, double tEnd)
{
    // In this method, we iterate over timesteps, doing the following for each:
    //   - update V using a forward Euler step
    //   - call ComputeExceptVoltage(t) to update the remaining state variables
    //     using backward Euler
    // Check length of time interval
    double _n_steps = (tEnd - tStart) / mDt;
    unsigned n_steps = (unsigned) round(_n_steps);
    assert(fabs(tStart+n_steps*mDt - tEnd) < 1e-12);
    
    // Initialise solution store
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(n_steps);
    solutions.rGetSolutions().push_back(rGetStateVariables());
    solutions.rGetTimes().push_back(tStart);
    
    // Loop over time
    double curr_time;
    for (unsigned i=0; i<n_steps; i++)
    {
        curr_time = tStart + i*mDt;
        
        // Compute next value of V
        UpdateTransmembranePotential(curr_time);
        
        // Compute other state variables
        ComputeOneStepExceptVoltage(curr_time);
        
        // Update solutions
        solutions.rGetSolutions().push_back(rGetStateVariables());
        solutions.rGetTimes().push_back(curr_time+mDt);
        
        // check gating variables are still in range
        VerifyGatingVariables();
    }
    
    return solutions;
}

template <unsigned SIZE>
void AbstractBackwardEulerCardiacCell<SIZE>::ComputeExceptVoltage(double tStart, double tEnd)
{
    // This method iterates over timesteps, calling ComputeExceptVoltage(t) at
    // each one, to update all state variables except for V, using backward Euler.
    // Check length of time interval
    double _n_steps = (tEnd - tStart) / mDt;
    unsigned n_steps = (unsigned) round(_n_steps);
    assert(fabs(tStart+n_steps*mDt - tEnd) < 1e-12);
    
    // Loop over time
    double curr_time;
    for (unsigned i=0; i<n_steps; i++)
    {
        curr_time = tStart + i*mDt;
        
        // Compute other state variables
        ComputeOneStepExceptVoltage(curr_time);
        
        // check gating variables are still in range
        VerifyGatingVariables();
    }
}


#endif /*ABSTRACTBACKWARDEULERCARDIACCELL_HPP_*/
