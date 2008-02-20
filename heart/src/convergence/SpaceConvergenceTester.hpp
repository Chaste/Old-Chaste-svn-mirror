#ifndef SPACECONVERGENCETESTER_HPP_
#define SPACECONVERGENCETESTER_HPP_


#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class SpaceConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM, PROBLEM_DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->MeshNum=0;
    }
    void UpdateConvergenceParameters()
    {
        this->MeshNum++;
    
    }
    bool GiveUpConvergence()
    {
        switch(DIM)
        {
            case 1:
            {
                return this->MeshNum>6;
                break;
            }
            case 2:
            {
                return this->MeshNum>6;
                break;
            }
            case 3:
            {
                return this->MeshNum>4;
                break;
            }
            default:
                assert(0);
                return true;//To keep Inntel compiler happy
        }
        return true;//To keep Intel compiler happy
    }
    double Abscissa()
    {
        unsigned mesh_size = (unsigned) pow(2, this->MeshNum+2); // number of elements in each dimension
        return this->mMeshWidth/(double) mesh_size;
    }
    
    int GetMeshNum()
    {
        return (int) this->MeshNum; //unsigned -> int is just cosmetic here.  (The test looks prettier).
    }
    double GetSpaceStep()
    {
        unsigned mesh_size = (unsigned) pow(2, this->MeshNum+2);// number of elements in each dimension
        double scaling = this->mMeshWidth/(double) mesh_size;
        return scaling;
    }
    
};

#endif /*SPACECONVERGENCETESTER_HPP_*/
