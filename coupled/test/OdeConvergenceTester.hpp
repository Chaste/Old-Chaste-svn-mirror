#ifndef ODECONVERGENCETESTER_HPP_
#define ODECONVERGENCETESTER_HPP_

#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class OdeConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->PdeTimeStep = 2e-2;
        this->OdeTimeStep = 1e-2;
    }
    void UpdateConvergenceParameters()
    {
        this->OdeTimeStep *= 0.5;
    
    }
    bool GiveUpConvergence()
    {
        return this->OdeTimeStep<=1e-8;
    }
    
    double Abscissa()
    {
        return this->OdeTimeStep;
    }
    
};
#endif /*ODECONVERGENCETESTER_HPP_*/
