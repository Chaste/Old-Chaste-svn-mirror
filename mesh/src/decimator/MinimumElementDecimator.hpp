#ifndef _MINIMUMELEMENTDECIMATOR_HPP_
#define _MINIMUMELEMENTDECIMATOR_HPP_


#include "Decimator.hpp"

///Note ELEMENT_DIM matches SPACE_DIM
template <unsigned SPACE_DIM>
class MinimumElementDecimator : public Decimator<SPACE_DIM>
{
protected:

    /// Returns the volume of the smallest element in the neighbourhood
    void CalculateLocalMeasure(NodeInfo<SPACE_DIM> *pNodeInfo, bool before)
    {
        double measure=INFINITY;
        Node<SPACE_DIM>* p_node=pNodeInfo->pGetNode();
        for (unsigned i=0; i<p_node->GetNumContainingElements();i++)
        {
        
            Element<SPACE_DIM,SPACE_DIM> *p_element=this->mpMesh->GetElement(p_node->GetNextContainingElementIndex());
            double volume=p_element->GetJacobianDeterminant();
            if (volume != 0.0 && volume < measure)
            {
                measure=volume;
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
            //A small element will be removed
            return this->mMeasureBefore;
        }
        //Otherwise a smaller element has been introduced
        return INFINITY;
    }
public:

};


#endif //_MINIMUMELEMENTDECIMATOR_HPP_
