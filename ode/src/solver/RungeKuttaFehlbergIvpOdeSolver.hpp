/**
 * Concrete RungeKutta2IvpOdeSolver class.
 */
#ifndef _RUNGEKUTTAFEHLBERGIVPODESOLVER_HPP_
#define _RUNGEKUTTAFEHLBERGIVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class RungeKuttaFehlbergIvpOdeSolver : public AbstractIvpOdeSolver
{
friend class TestRungeKuttaFehlbergIvpOdeSolver;    
    
private:
    /*
     * All these are here for more efficient memory allocation, rather than 
     * because they need to be member variables...
     */
    double m1932o2197;
    double m7200o2197;
    double m7296o2197;
    double m12o13;
    double m439o216;
    double m3680o513;
    double m845o4104;
    double m8o27;
    double m3544o2565;
    double m1859o4104;
    double m1o360;
    double m128o4275;
    double m2197o75240;
    double m2o55;
    double m25o216;
    double m1408o2565;
    double m2197o4104;
    std::vector<double> mError;
    std::vector<double> mk1;
    std::vector<double> mk2;
    std::vector<double> mk3;
    std::vector<double> mk4;
    std::vector<double> mk5;
    std::vector<double> mk6;
    std::vector<double> myk2;
    std::vector<double> myk3;
    std::vector<double> myk4;
    std::vector<double> myk5;
    std::vector<double> myk6;
    
protected:
    /**
     * Method that actually performs the solving on behalf of the public Solve methods.
     * 
     * @param rSolution  an ode solution to input data into if requited
     * @param pAbstractOdeSystem  the system to solve
     * @param rCurrentYValues  the current (initial) state; results will also be returned
     *     in here
     * @param rWorkingMemory  working memory; same size as rCurrentYValues
     * @param startTime  initial time
     * @param endTime  time to solve to
     * @param timeStep  dt
     * @param outputSolution whether to output into rSolution (or save time by not doing)
     */
    void InternalSolve(OdeSolution& rSolution,
                               AbstractOdeSystem* pAbstractOdeSystem,
                               std::vector<double>& rCurrentYValues,
                               std::vector<double>& rWorkingMemory,
                               double startTime,
                               double endTime,
                               double maxTimeStep,
                               double minTimeStep,
                               double tolerance,
                               bool outputSolution);
                               
    /**
     * Calculate the next time step. Updates the mError vector with current error.
     */
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                     double timeStep,
                                     double time,
                                     std::vector<double>& currentYValues,
                                     std::vector<double>& nextYValues);
    
    /**
     * @param rCurrentStepSize The current step size being used (returns answer via this reference)
     * @param rError The error in the approximation at this time step
     * @param rTolerance The tolerance required
     * @param rMaxTimeStep The maximum timestep to be used
     * @param rMinTimeStep The minimum timestep to be used (to prevent huge loops)
     */                                 
    void AdjustStepSize(double& rCurrentStepSize, const double& rError, const double& rTolerance, 
                                const double& rMaxTimeStep, const double& rMinTimeStep);                               
                                     
public:

    OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem,
                              std::vector<double>& rYValues,
                              double startTime,
                              double endTime,
                              double timeStep,
                              double ignoredSamplingTime);
                              
    void Solve(AbstractOdeSystem* pAbstractOdeSystem,
                       std::vector<double>& rYValues,
                       double startTime,
                       double endTime,
                       double timeStep);                       
                       
    RungeKuttaFehlbergIvpOdeSolver()
      : m1932o2197(1932.0/2197.0),// you need the .0 s - caused me no end of trouble.
        m7200o2197(7200.0/2197.0),
        m7296o2197(7296.0/2197.0),
        m12o13(12.0/13.0),
        m439o216(439.0/216.0),
        m3680o513(3680.0/513.0),
        m845o4104(845.0/4104.0),
        m8o27(8.0/27.0),
        m3544o2565(3544.0/2565.0),
        m1859o4104(1859.0/4104.0),
        m1o360(1.0/360.0),
        m128o4275(128.0/4275.0),
        m2197o75240(2197.0/75240.0),
        m2o55(2.0/55.0),
        m25o216(25.0/216.0),
        m1408o2565(1408.0/2565.0),
        m2197o4104(2197.0/4104.0)
    {
    };
};

#endif //_RUNGEKUTTAFEHLBERGIVPODESOLVER_HPP_
