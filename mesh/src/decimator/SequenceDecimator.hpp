#ifndef _SEQUENCE_DECIMATOR_HPP_
#define _SEQUENCE_DECIMATOR_HPP_


#include "Decimator.hpp"

///Note ELEMENT_DIM matches SPACE_DIM
template <unsigned SPACE_DIM>
class SequenceDecimator : public Decimator<SPACE_DIM>
{
protected:

    void CalculateLocalMeasure(NodeInfo<SPACE_DIM> *pNodeInfo, bool before)
    {
    
    
        double measure=pNodeInfo->GetIndex();
        if (before)
        {
            this->mMeasureBefore=measure;
        }
        else
        {
            this->mMeasureAfter=measure;
        }
        
    }
    
    //double CalculateScore returns mMeasureAfter as in the base class
public:

};


#endif //_SEQUENCE_DECIMATOR_HPP_
