#ifndef _ABSTRACTODESYSTEMWITHANALYTICJACOBIAN_HPP_
#define _ABSTRACTODESYSTEMWITHANALYTICJACOBIAN_HPP_

#include <vector>

#include "Exception.hpp"
#include "AbstractOdeSystem.hpp"

/**
 * Abstract Analytic Jacobian
 *
 * Represents an ODE system with an analytic Jacobian available,
 * which can be computed using the method AnalyticJacobian.
 */
class AbstractOdeSystemWithAnalyticJacobian : public AbstractOdeSystem
{
protected:

public:
    AbstractOdeSystemWithAnalyticJacobian(unsigned numberOfStateVariables = 0)
        : AbstractOdeSystem(numberOfStateVariables)
    {
        mUseAnalytic = true;
    }
    
    /**
     * Compute the analytic Jacobian matrix of the ODE system.
     *
     * @param rSolutionGuess  the current guess at the solution for this time step
     * @param jacobian  will be filled in with the Jacobian matrix entries
     * @param time  the current simulation time
     * @param timeStep  the time step in use by the integrator at present
     */
    virtual void AnalyticJacobian(const std::vector<double>& rSolutionGuess,
                                  double** jacobian,
                                  double time,
                                  double timeStep) = 0;
    
};

#endif //_ABSTRACTODESYSTEMWITHANALYTICJACOBIAN_HPP_
