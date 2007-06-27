#ifndef _QUALITYDECIMATOR_HPP_
#define _QUALITYDECIMATOR_HPP_


#include "Decimator.hpp"

///Note ELEMENT_DIM matches SPACE_DIM
template <unsigned SPACE_DIM>
class QualityDecimator : public Decimator<SPACE_DIM>
{
protected:

    void CalculateLocalMeasure(NodeInfo<SPACE_DIM> *pNodeInfo, bool before)
    {
        double measure=INFINITY;
        Node<SPACE_DIM> *p_node=pNodeInfo->pGetNode();
        for (typename Node<SPACE_DIM>::ContainingElementIterator it = p_node->ContainingElementsBegin();
             it != p_node->ContainingElementsEnd();
             ++it)
        {
            Element<SPACE_DIM,SPACE_DIM> *p_element=this->mpMesh->GetElement(*it);
            double volume=p_element->GetJacobianDeterminant();
            if (volume != 0.0)
            {
                double quality=p_element->CalculateQuality();
                if (quality < measure)
                {
                    measure=quality;
                }
            }
        }
        if (before)
        {
            this->mMeasureBefore=measure;
        }
        else
        {
            this->mMeasureAfter=measure;
        }
    }
    
    double CalculateScore()
    {
        if (this->mMeasureBefore < this->mMeasureAfter)
        {
            //A bad quality element will be removed
            return this->mMeasureBefore;
        }
        //Otherwise a worse element has been introduced
        return INFINITY;
    }
public:

};


#endif //_QUALITYDECIMATOR_HPP_
