#ifndef ABSTRACTBACKWARDEULERCARDIACCELL_HPP_
#define ABSTRACTBACKWARDEULERCARDIACCELL_HPP_

#include "AbstractCardiacCell.hpp"
//#include "CardiacNewtonSolver.hpp"

#include <cassert>

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
                        unsigned int numberOfStateVariables,
                        unsigned int voltageIndex,
                        double dt,
                        AbstractStimulusFunction* intracellularStimulus,
                        AbstractStimulusFunction* extracellularStimulus = NULL)
            : AbstractCardiacCell(NULL,
                                  numberOfStateVariables,
                                  voltageIndex,
                                  dt,
                                  intracellularStimulus,
                                  extracellularStimulus)
    {
    }
    
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
     * We override the default Compute* method to make them pure virtual, since the
     * default implementations are not suitable - we don't have access to a standard
     * ODE solver.
     */    
    virtual OdeSolution Compute(double tStart, double tEnd)=0;
    virtual OdeSolution ComputeExceptVoltage(double tStart, double tEnd)=0;
    
#define COVERAGE_IGNORE
    /**
     * This function should never be called - the cell class incorporates its own solver.
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
    {
        assert(0);
    }
#undef COVERAGE_IGNORE
};

#endif /*ABSTRACTBACKWARDEULERCARDIACCELL_HPP_*/
