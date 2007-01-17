#ifndef _EFFICIENTEULERIVPODESOLVER_HPP_
#define _EFFICIENTEULERIVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>
#include <cassert>


/**
 * A version of the EulerIvpOdeSolver that makes more efficient use of
 * memory, by not repeatedly creating vectors unnecessarily.
 *
 * The aim is to see if the use of lookup tables can lead to a greater
 * speedup than seen thus far.
 *
 * Solves a system of ODEs using the Forward Euler method.
 *
 * To be used in the form:
 *
 * EfficientEulerIvpOdeSolver my_solver;
 *
 * OdeSolution solution = my_solver.Solve(pMyOdeSystem, yInit,
 *                                        StartTime, EndTime, TimeStep, SamplingTime);
 *
 * See also documentation for EfficientEulerIvpOdeSolver::Solve()
 *
 */
class EfficientEulerIvpOdeSolver : public AbstractOneStepIvpOdeSolver
{
public:
    EfficientEulerIvpOdeSolver()
    {} //Constructor-does nothing
    
    virtual OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem,
                  std::vector<double>& rInitialYValues,
                  double startTime,
                  double endTime,
                  double timeStep,
                  double timeSampling);

    virtual void Solve(AbstractOdeSystem* pAbstractOdeSystem,
                       std::vector<double>& rInitialYValues,
                       double startTime,
                       double endTime,
                       double timeStep);                  
    
    virtual void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                     double timeStep,
                     double time,
                     std::vector<double>& rCurrentYValues,
                     std::vector<double>& rNextYValues);

    /** 
     * We need this method since it's a pure virtual method in the base
     * class.  But don't use it.
     */
    virtual std::vector<double> CalculateNextYValue(
        AbstractOdeSystem* pAbstractOdeSystem,
        double timeStep,
        double time,
        std::vector<double> rCurrentYValues)
    {
        assert(0);
        return rCurrentYValues;
    }
    
    virtual ~EfficientEulerIvpOdeSolver()
    {}
    
};

#endif //_EFFICIENTEULERIVPODESOLVER_HPP_
