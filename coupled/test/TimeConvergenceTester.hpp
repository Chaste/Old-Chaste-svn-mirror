#ifndef TIMECONVERGENCETESTER_HPP_
#define TIMECONVERGENCETESTER_HPP_

#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class TimeConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->mPdeTimeStep = 0.04;
        this->mOdeTimeStep = this->mPdeTimeStep;
    }
    void UpdateConvergenceParameters()
    {
        std::cout << "Yo" << "\n";
        this->mPdeTimeStep *= 0.5;
        this->mOdeTimeStep = this->mPdeTimeStep;
    
    }
    bool GiveUpConvergence()
    {
        return this->mPdeTimeStep<=1e-8;
    }
    
};
#endif /*TIMECONVERGENCETESTER_HPP_*/
