#ifndef PDECONVERGENCETESTER_HPP_
#define PDECONVERGENCETESTER_HPP_

#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class PdeConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->PdeTimeStep = 0.04;
        this->OdeTimeStep = 0.0025;
    }
    void UpdateConvergenceParameters()
    {
        this->PdeTimeStep *= 0.5;
    }
    bool GiveUpConvergence()
    {
        return this->PdeTimeStep<=1e-8;
    }
    double Abscissa()
    {
        return this->PdeTimeStep;
    }
};
#endif /*PDECONVERGENCETESTER_HPP_*/
