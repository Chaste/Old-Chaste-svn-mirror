#ifndef SPACECONVERGENCETESTER_HPP_
#define SPACESCONVERGENCETESTER_HPP_


#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class SpaceConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->mMeshNum=0;
    }
    void UpdateConvergenceParameters()
    {
        this->mMeshNum++;
    
    }
    bool GiveUpConvergence()
    {
        switch(DIM)
        {
            case 1:
            {
                return this->mMeshNum>6;
                break;
            }
            case 2:
            {
                return this->mMeshNum>4;
                break;
            }
            case 3:
            {
                assert(0);
                break;
            }
            default:
                assert(0);
        }
    }
    double Abscissa()
    {
        unsigned mesh_size = (unsigned) pow(2, this->mMeshNum+2); // number of elements in each dimension
        return mesh_width/(double) mesh_size;
    }
};

#endif /*SPACECONVERGENCETESTER_HPP_*/
