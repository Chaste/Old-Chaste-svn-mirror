#ifndef SPACECONVERGENCETESTER_HPP_
#define SPACECONVERGENCETESTER_HPP_


#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class SpaceConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
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
                return this->MeshNum>1;
                break;
            }
            case 2:
            {
                return this->MeshNum>5;
                break;
            }
            case 3:
            {
                assert(0);
                return true;//To keep Intel compiler happy
                break;
            }
            default:
                assert(0);
                return true;//To keep Intel compiler happy
        }
        return true;//To keep Intel compiler happy
    }
    double Abscissa()
    {
        unsigned mesh_size = (unsigned) pow(2, this->MeshNum+2); // number of elements in each dimension
        return mesh_width/(double) mesh_size;
    }
};

#endif /*SPACECONVERGENCETESTER_HPP_*/
