#ifndef SPACECONVERGENCETESTER_HPP_
#define SPACESCONVERGENCETESTER_HPP_


#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class SpaceConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        std::cout << "Yep";
        this->mMeshNum=1;
    }
    void UpdateConvergenceParameters()
    {
        this->mMeshNum++;
    
    }
    bool GiveUpConvergence()
    {
        return this->mMeshNum>6;
    }
};

#endif /*SPACECONVERGENCETESTER_HPP_*/
