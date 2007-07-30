#ifndef _VECTOR_FUNCTION_DECIMATOR_HPP_
#define _VECTOR_FUNCTION_DECIMATOR_HPP_


#include "Decimator.hpp"
template <unsigned SPACE_DIM>
class VectorFunctionNodeInfo : public NodeInfo<SPACE_DIM>
{
public:
    VectorFunctionNodeInfo(Node<SPACE_DIM> *pNode, unsigned position)
    : NodeInfo<SPACE_DIM>(pNode, position)
    {
         
    }

private:
    bool
    TieBreakGreaterThan(const NodeInfo<SPACE_DIM>* other) const
    {
        return  (this->mNeighbourhoodVolume < other->GetNeighbourhoodVolume());
        
    }
};
///Note ELEMENT_DIM matches SPACE_DIM
template <unsigned SPACE_DIM>
class VectorFunctionDecimator : public Decimator<SPACE_DIM>
{
private:
    std::vector<c_vector<double, SPACE_DIM> > mPayload;
    out_stream mFibreFile;
    c_vector<double, SPACE_DIM> mVectorMeasureBefore, mVectorMeasureAfter;
    
    NodeInfo<SPACE_DIM> *CreateNodeInfo(Node<SPACE_DIM> *p_node, unsigned position)
    {
        return new VectorFunctionNodeInfo<SPACE_DIM>(p_node, position);
    }
protected:
    void CalculateLocalMeasure(NodeInfo<SPACE_DIM> *pNodeInfo, bool before)
    {
        c_vector<double, SPACE_DIM> measure=zero_vector<double>(SPACE_DIM);
        Node<SPACE_DIM> *p_node=pNodeInfo->pGetNode();
        for (typename Node<SPACE_DIM>::ContainingElementIterator it = p_node->ContainingElementsBegin();
             it != p_node->ContainingElementsEnd();
             ++it)
        {
            Element<SPACE_DIM,SPACE_DIM> *p_element = this->mpMesh->GetElement(*it);
            double volume = p_element->GetJacobianDeterminant();
            if (volume != 0.0)
            {
                c_vector<double, SPACE_DIM> weight=zero_vector<double>(SPACE_DIM);
                for (unsigned j=0; j<=SPACE_DIM;j++)
                {
                    unsigned index=p_element->GetNodeGlobalIndex(j);
                    if (before == false && index==(unsigned)p_node->GetIndex())
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
            this->mVectorMeasureBefore=measure;
        }
        else
        {
            this->mVectorMeasureAfter=measure;
        }
    }
    
    double CalculateScore()
    {
    
    
        double measure=norm_2(mVectorMeasureBefore-mVectorMeasureAfter);
        //Is the measure negligible when compared to the volume?
        if (measure/this->mNeighbourhoodVolume < 1e-3)
        {
            //std::cout<<measure<<"\t"<<this->mNeighbourhoodVolume<<"\t"<<"Small\n";
            //measure += this->mNeighbourhoodVolume;
            //std::cout<<measure<<"\t"<<this->mNeighbourhoodVolume<<"\t"<<"Small after fudge\n";
            
        }
        else
        {
            //std::cout<<measure<<"\t"<<this->mNeighbourhoodVolume<<"\t"<<"Large\n";
        }
        return measure;
    }
    
    virtual void OpenAnimationFiles(std::string filePathName)
    {
        OutputFileHandler handler(filePathName+"/vis_results/",true);
        
        this->mNodeFile=handler.OpenOutputFile("results.viznodes");
        this->mElementFile=handler.OpenOutputFile("results.vizelements");
        mFibreFile=handler.OpenOutputFile("results.vizfibres");
    }
    
    virtual void CloseAnimationFiles()
    {
        this->mNodeFile->close();
        this->mElementFile->close();
        mFibreFile->close();
        
    }
    
    virtual void WriteVisualiseFiles(double time)
    {
        Decimator<SPACE_DIM>::WriteVisualiseFiles(time);
        
        (*mFibreFile)<<time<<"\t";
        
        for (unsigned i=0; i<(unsigned)this->mpMesh->GetNumAllNodes();i++)
        {
            Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(i);
            
            if (p_node->IsDeleted() == false)
            {
                (*mFibreFile) << mPayload[i](0) <<"\t"<<mPayload[i](1)<<"\t";
            }
        }
        (*mFibreFile) <<"\n";
    }
    
public:
    void Initialise(ConformingTetrahedralMesh<SPACE_DIM, SPACE_DIM> *pMesh, std::vector<c_vector<double,SPACE_DIM> > payload)
    {
        mPayload=payload;
        Decimator<SPACE_DIM>::Initialise(pMesh);
    }
};


#endif //_VECTOR_FUNCTION_DECIMATOR_HPP_
