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

public:
    
    BackwardEulerIvpOdeSolver()
    {
        // default epsilon
       mEpsilon = 1e-6;
    }
    
    void SetEpsilonForNumericalJacobian(double epsilon)
    {
        assert(epsilon > 0);
        mEpsilon = epsilon;
    }
    
    double GetMEplison()
    {
        return mEpsilon;
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
