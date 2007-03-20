#ifndef _RANDOM_DECIMATOR_HPP_
#define _RANDOM_DECIMATOR_HPP_


#include "Decimator.hpp"
#include "RandomNumberGenerator.hpp"
///Note ELEMENT_DIM matches SPACE_DIM
template <unsigned SPACE_DIM>
class RandomDecimator : public Decimator<SPACE_DIM>
{
protected:

    void CalculateLocalMeasure(NodeInfo<SPACE_DIM> *pNodeInfo, bool before)
    {
        if (before)
        {
            this->mMeasureBefore=RandomNumberGenerator::Instance()->ranf();
        }
        else
        {
            this->mMeasureAfter=this->mMeasureBefore;
        }
        
    }
    
    //double CalculateScore returns mMeasureAfter as in the base class
    
};


#endif //_RANDOM_DECIMATOR_HPP_
