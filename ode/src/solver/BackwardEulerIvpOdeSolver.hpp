/**
 * Concrete BackwardEulerIvpOdeSolver class.
 */
#ifndef BACKWARDEULERIVPODESOLVER_HPP_
#define BACKWARDEULERIVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"
#include <cassert>
#include <vector>

class BackwardEulerIvpOdeSolver  : public AbstractOneStepIvpOdeSolver
{
private:
    /** the epsilon to use in calculating numerical jacobians */
    double mEpsilon;
    bool mForceUseOfNumericalJacobian;

public:
    
    BackwardEulerIvpOdeSolver()
    {
        // default epsilon
        mEpsilon = 1e-6;
        mForceUseOfNumericalJacobian = false;
    }
    
    void SetEpsilonForNumericalJacobian(double epsilon)
    {
        assert(epsilon > 0);
        mEpsilon = epsilon;
    }
     
    /** Force the solver to use the numerical Jacobian even if the 
     *  ode system is one with an analytical jacobian provided
     */
    void ForceUseOfNumericalJacobian()
    {
        mForceUseOfNumericalJacobian = true;
    }
    
    virtual ~BackwardEulerIvpOdeSolver()
    {
    }   
    
    
    
    std::vector<double> CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                            double timeStep,
                                            double time,
                                            std::vector<double> currentYValue);                                
                                            
};


#endif /*BACKWARDEULERIVPODESOLVER_HPP_*/
