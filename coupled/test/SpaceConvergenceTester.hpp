#ifndef SPACECONVERGENCETESTER_HPP_
#define SPACESCONVERGENCETESTER_HPP_


#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class SpaceConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
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
    double Abscissa()
    {
        unsigned mesh_size = (unsigned) pow(2, this->mMeshNum+2); // number of elements in each dimension
        return mesh_width/(double) mesh_size;
    }
};

#endif /*SPACECONVERGENCETESTER_HPP_*/
