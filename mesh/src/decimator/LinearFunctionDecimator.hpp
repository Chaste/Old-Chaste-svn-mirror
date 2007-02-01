#ifndef _LINEAR_FUNCTION_DECIMATOR_HPP_
#define _LINEAR_FUNCTION_DECIMATOR_HPP_


#include "Decimator.hpp"

///Note ELEMENT_DIM matches SPACE_DIM
template <int SPACE_DIM>
class LinearFunctionDecimator : public Decimator<SPACE_DIM>
{
private:
    std::vector<double> mPayload; 
protected:
    
    void CalculateLocalMeasure(NodeInfo<SPACE_DIM> *pNodeInfo, bool before)
    {
        double measure=0.0;
        Node<SPACE_DIM> *p_node=pNodeInfo->pGetNode();
        for (unsigned i=0; i<p_node->GetNumContainingElements();i++)
        {
            Element<SPACE_DIM,SPACE_DIM> *p_element=
                    this->mpMesh->GetElement(p_node->GetNextContainingElementIndex());
            double volume = p_element->GetJacobianDeterminant();
            if (volume != 0.0)
            {
                //Need to watch out for where it's going
                double weight=0.0;
                for (int j=0; j<=SPACE_DIM;j++)
                {
                    int index=p_element->GetNodeGlobalIndex(j);
                    if (before == false && index==p_node->GetIndex())
                    {
                        index=pNodeInfo->GetPossibleTargetIndex();
                    }
                    weight+=mPayload[index];
                }
                measure += weight*volume;
            }
        }
        if (before)
        {
            this->mMeasureBefore=measure;
        } else {
            this->mMeasureAfter=measure;
        }
    }
    
    double CalculateScore()
    {
        //std::cout << "Diff = "<<fabs(this->mMeasureBefore-this->mMeasureAfter);
        //std::cout << "\t vol= "<<this->mNeighbourhoodVolume<<"\n";
        
        double measure=fabs(this->mMeasureBefore-this->mMeasureAfter);
        
        //Is the measure negligible when compared to the volume?
        if (measure/this->mNeighbourhoodVolume < 1e-3){
            //std::cout<<"NEG\n";
            measure += this->mNeighbourhoodVolume*1e-6;
        
        }
        return fabs(measure);
    }
    
    int GetTag(int index)
    {
        if (this->mpMesh->GetNode(index)->IsBoundaryNode())
        {
            return 0;
        }
        if (mPayload[index]<0)
        {   
            return 1;
        }
        return 2; 
    }
public:
    void Initialise(ConformingTetrahedralMesh<SPACE_DIM, SPACE_DIM> *pMesh, std::vector<double> payload)
    {
        mPayload=payload;
        Decimator<SPACE_DIM>::Initialise(pMesh);
    }
};


#endif //_LINEAR_FUNCTION_DECIMATOR_HPP_
