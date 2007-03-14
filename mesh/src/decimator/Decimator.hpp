#ifndef _DECIMATOR_HPP_
#define _DECIMATOR_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "NodeInfo.hpp"
#include "Exception.hpp"
#include <vector>
#include <algorithm>

/*** \todo
* makeheap efficiency
* inner product interpolation
* zero area surface?
*/

///Note ELEMENT_DIM matches SPACE_DIM
template <unsigned SPACE_DIM>
class Decimator
{
protected:
    ConformingTetrahedralMesh<SPACE_DIM, SPACE_DIM> *mpMesh;
    std::vector<NodeInfo<SPACE_DIM> *> mQueue;
    double mThresholdScore;
    double mVolumeLeakage;
    double mMeasureBefore;
    double mMeasureAfter;
    double mNeighbourhoodVolume;
    out_stream mNodeFile, mElementFile, mFibreFile;
    
    void Initialise(NodeInfo<SPACE_DIM> *pNodeInfo)
    {
        //Pre-requesite mQueue is such that mQueue[i]->GetIndex() == i
        Node<SPACE_DIM> *pNode=pNodeInfo->mpNode;
        for (unsigned i=0; i<pNode->GetNumContainingElements(); i++)
        {
            unsigned element_index=pNode->GetNextContainingElementIndex();
            Element <SPACE_DIM, SPACE_DIM> *p_element= mpMesh->GetElement(element_index);
            for (unsigned j=0; j<SPACE_DIM+1; j++)
            {
                unsigned index=p_element->GetNodeGlobalIndex(j);
                pNodeInfo->AddToNeighbourNodes(mQueue[index]);
            }
        }
       
    }
    
    virtual void CalculateLocalMeasure(NodeInfo<SPACE_DIM> *pNodeInfo, bool before)
    {
        double measure=0.0;
        if (before) measure=0.0;
        Node<SPACE_DIM> *p_node=pNodeInfo->mpNode;
        for (unsigned i=0; i<p_node->GetNumContainingElements();i++)
        {
            Element<SPACE_DIM,SPACE_DIM> *p_element=mpMesh->GetElement(p_node->GetNextContainingElementIndex());
            measure += p_element->GetJacobianDeterminant();
        }
        if (before)
        {
            mMeasureBefore=measure;
        } else {
            mMeasureAfter=measure;
        }
    }
    
    double CalculateNeighbourhoodVolume(Node<SPACE_DIM> *pNode)
    {
        double measure=0.0;
        for (unsigned i=0; i<pNode->GetNumContainingElements();i++)
        {
            Element<SPACE_DIM,SPACE_DIM> *p_element=mpMesh->GetElement(pNode->GetNextContainingElementIndex());
            measure += p_element->GetJacobianDeterminant();
        }
  
        return measure;
    }
    
    virtual double CalculateScore()
    {
        //Double check for volume leakage in base decimator
        assert( fabs(mMeasureBefore-mMeasureAfter)/fabs(mMeasureBefore+mMeasureAfter) 
                <  mVolumeLeakage);

        return mMeasureBefore;
    }
    
    double CalculateBoundaryScore(double neighbourhoodBefore, double neighbourhoodAfter)
    {
        double change=fabs(neighbourhoodAfter-neighbourhoodBefore);
        change /= (neighbourhoodBefore+neighbourhoodAfter);
        if (change<mVolumeLeakage)
        {
            return 0.0;
        }
        return INFINITY;
    }
    void Rescore(NodeInfo<SPACE_DIM> *pNodeInfo)
    {
        pNodeInfo->mScore=INFINITY;
        //pNodeInfo->mPossibleTargetIndex=pNodeInfo->GetIndex();
        CalculateLocalMeasure(pNodeInfo, true);
        Point<SPACE_DIM> point=pNodeInfo->mpNode->GetPoint();


        mNeighbourhoodVolume=CalculateNeighbourhoodVolume(pNodeInfo->mpNode);
        

        for (unsigned i=0; i<pNodeInfo->GetNumNeighbourNodes();i++)
        {
            NodeInfo<SPACE_DIM> *p_neighbour=pNodeInfo->GetNextNeighbourNode();
            pNodeInfo->mPossibleTargetIndex=p_neighbour->GetIndex();
            try
            {
                mpMesh->MoveMergeNode(pNodeInfo->GetIndex(),
                                pNodeInfo->mPossibleTargetIndex, false);
        		double score=0.0;
                if (pNodeInfo->mpNode->IsBoundaryNode())
                {
                    double neighbourhood_volume
                                =CalculateNeighbourhoodVolume(pNodeInfo->mpNode);
                    score += CalculateBoundaryScore(mNeighbourhoodVolume,
                            neighbourhood_volume);
                }
                if (score != INFINITY)
                {
                    CalculateLocalMeasure(pNodeInfo, false);
                    score+=CalculateScore();
                }
                if (score < pNodeInfo->mScore)
                {
                    pNodeInfo->mScore=score;
                    pNodeInfo->mpBestNeighbourNode=p_neighbour;
                }
                if (score != INFINITY && score == pNodeInfo->mScore)
                {
                    //Always select predictable index for merger
                    if  (pNodeInfo->mpBestNeighbourNode->mpNode->GetIndex() < p_neighbour->mpNode->GetIndex())
                    {
                        pNodeInfo->mpBestNeighbourNode=p_neighbour;
                    }
                }
                
            }
            catch (Exception e)
            {
                //If the move is not feasible then we ignore it
  
            }
        }
        mpMesh->SetNode(pNodeInfo->GetIndex(), point);
        
    }
    
    
    
    void ActivateOnce()
    {
        NodeInfo<SPACE_DIM> *p_moving_node_info=mQueue[0];
        unsigned moving_node_index=p_moving_node_info->GetIndex();
        
        NodeInfo<SPACE_DIM> *p_target_node_info=p_moving_node_info->mpBestNeighbourNode;
        unsigned target_node_index=p_target_node_info->GetIndex();
        
        
        
        //std::cout<<"Moving "<<moving_node_index<<" to "<<target_node_index<<"\n";
        mpMesh->MoveMergeNode(moving_node_index, target_node_index);
        
        //Update the neighbour information
        //The target node will have all the moving nodes neighbours added to its set
        //All other neighbours of the moving node will have the target node added and the moving
        //node removed.
        
        //Add neighbour nodes to target
        for (unsigned i=0; i<p_moving_node_info->GetNumNeighbourNodes(); i++)
        {
            p_target_node_info->AddToNeighbourNodes(p_moving_node_info->GetNextNeighbourNode());
        }
        
        /*
        //Add neighbour elements to target
        //This is unnecessary since the moving node has already de-registered from
        // its containing elements (when marked for deletion  
        Node<SPACE_DIM> *p_moving_node=p_moving_node_info->mpNode;
        Node<SPACE_DIM> *p_target_node=p_target_node_info->mpNode;
        for (unsigned i=0; i<p_moving_node->GetNumContainingElements(); i++)
        {
            p_target_node->AddElement(                 p_moving_node->GetNextContainingElementIndex());
        }
        */
        //Swap moving node for target in all neigbours
        for (unsigned i=0; i<p_moving_node_info->GetNumNeighbourNodes(); i++)
        {
            NodeInfo<SPACE_DIM> *p_neighbour_info=p_moving_node_info->GetNextNeighbourNode();
            p_neighbour_info->RemoveFromNeighbourNodes(p_moving_node_info);
            p_neighbour_info->AddToNeighbourNodes(p_target_node_info);
        }
        
        
        
        
        for (unsigned i=0; i<p_moving_node_info->GetNumNeighbourNodes(); i++)
        {
            //Find the neighbour's information
            NodeInfo<SPACE_DIM> *p_neighbour_info=p_moving_node_info->GetNextNeighbourNode();
            Rescore(p_neighbour_info);
        }
        //Remove the moving node information from further consideration
        std::pop_heap(mQueue.begin(), mQueue.end(), CompNodeInfo<SPACE_DIM>());
        mQueue.pop_back();
        delete(p_moving_node_info);
        
        //This is potentially inefficient since we've only changed the odd node here and there
        std::make_heap(mQueue.begin(),mQueue.end(),CompNodeInfo<SPACE_DIM>());
        
    }
    
    
    
    virtual void OpenAnimationFiles(std::string filePathName)
    {
        OutputFileHandler handler("");
        mNodeFile=handler.OpenOutputFile(filePathName+".viznodes");
        mElementFile=handler.OpenOutputFile(filePathName+".vizelements");
    }
    
    virtual void CloseAnimationFiles()
    {
        mNodeFile->close();
        mElementFile->close();
    }
   
    
    virtual void WriteVisualiseFiles(double time)
    {
        (*mNodeFile)<<time<<"\t";
        if (SPACE_DIM==2)
        {
            (*mElementFile)<<time<<"\t";
        }
        NodeMap node_map(mpMesh->GetNumAllNodes());
        
        unsigned new_index=0;
        for (unsigned i=0; i<(unsigned)mpMesh->GetNumAllNodes();i++)
        {
            Node<SPACE_DIM>* p_node = mpMesh->GetNode(i);
        
            if (p_node->IsDeleted() == false)
            {
                for (unsigned j=0; j<SPACE_DIM; j++)
                {
                    (*mNodeFile)<<p_node->GetPoint()[j]<<"\t";
                }
                node_map.SetNewIndex(i,new_index++);
                if (SPACE_DIM == 2)
                {
                    unsigned tag=GetTag(i);
                    (*mNodeFile)<<tag<<"\t";
                }
            }
            else
            {
                node_map.SetDeleted(i);
            }
        }
        (*mNodeFile)<<"\n";
        if (SPACE_DIM==2)
        {
            for (unsigned i=0; i<(unsigned)mpMesh->GetNumAllElements();i++)
            {
                Element<SPACE_DIM,SPACE_DIM> *element=mpMesh->GetElement(i);
                if (element->IsDeleted() == false)
                {
                    for (unsigned j=0; j<SPACE_DIM+1; j++)
                    {
                        unsigned old_index=element->GetNodeGlobalIndex(j);
                        (*mElementFile)<<node_map.GetNewIndex(old_index)<<"\t";
                    }
                }
            }
            (*mElementFile)<<"\n";
        }
        
    }
       
    virtual unsigned GetTag(unsigned index)
    {
        if (mpMesh->GetNode(index)->IsBoundaryNode())
        {
            return 0;
        }
        //else
        return 1; 
    }
public:
    void SetThreshold(double threshold)
    {
        mThresholdScore = threshold;
    }
    double GetThreshold()
    {
        return mThresholdScore;
    }
    
    void SetVolumeLeakage(double leakage)
    {
        mVolumeLeakage= leakage;
    }
    double GetVolumeLeakage()
    {
        return mVolumeLeakage;
    }
    
    void Decimate()
    {
        while (mQueue[0]->mScore < mThresholdScore)
        {
            ActivateOnce();
            //Interrogate();
        }
    }
    void DecimateAnimate(std::string filePathName, unsigned step=1)
    {
        unsigned time=0;
        assert(SPACE_DIM<=2);
        OpenAnimationFiles(filePathName);
        
        WriteVisualiseFiles(time);
        while (mQueue[0]->mScore < mThresholdScore)
        {
            ActivateOnce();
            time++;
            if ( time%step == 0){
                WriteVisualiseFiles(time);
            }
        }
        
        //Make sure that the final step is always shown
        if (step != 1 && time%step != 0){
            WriteVisualiseFiles(time);
        }
        CloseAnimationFiles();
        
    }
    
    
    void Initialise(ConformingTetrahedralMesh<SPACE_DIM, SPACE_DIM> *pMesh)
    {
        mpMesh = pMesh;
        //mVolume= pMesh->CalculateMeshVolume();
        mThresholdScore = INFINITY;
        mVolumeLeakage = 1e-5;
        mQueue.reserve(mpMesh->GetNumAllNodes());
        for (unsigned i=0; i<(unsigned)mpMesh->GetNumAllNodes();i++)
        {
            NodeInfo<SPACE_DIM> *p_node_info=new NodeInfo<SPACE_DIM>(mpMesh->GetNode(i));
            mQueue.push_back(p_node_info);
        }
        for (unsigned i=0; i<(unsigned)mpMesh->GetNumAllNodes();i++)
        {
            Initialise(mQueue[i]);
        }
        Rescore();
    }
    void Rescore()
    {
        for (unsigned i=0; i<(unsigned)mpMesh->GetNumAllNodes();i++)
        {
            Rescore(mQueue[i]);
        }
        std::make_heap(mQueue.begin(),mQueue.end(),CompNodeInfo<SPACE_DIM>());
    }
    virtual ~Decimator()
    {
        for (unsigned i=0; i<mQueue.size(); i++)
        {
            delete mQueue[i];
        }
    }
    
    void Interrogate()
    {
        std::cout<<"Number of nodes in queue is now "<<mQueue.size()<<"\t";
        std::cout<<"Lowest score is "<<mQueue[0]->mScore<<"\n";
        //return;
        std::vector<NodeInfo<SPACE_DIM> *> queue=mQueue;
        while ( !queue.empty() )
        {
            NodeInfo<SPACE_DIM> *p_node_info=queue[0];
            std::cout << "Score:"<<p_node_info->GetScore()<<"\t";
            std::cout << "Index:"<<p_node_info->GetIndex()<<"\t";
            std::cout << "Best:";
            if (p_node_info->mpBestNeighbourNode != 0)
            {
                std::cout<<p_node_info->mpBestNeighbourNode->GetIndex()<<"\n";
            }
            else
            {
                std::cout<<"(not set)\n";
            }
            for (unsigned i=0;i<p_node_info->GetNumNeighbourNodes();i++)
            {
                std::cout<<" "<<p_node_info->GetNextNeighbourNode()->GetIndex();
            }
            std::cout<<"\n";
            std::pop_heap(queue.begin(), queue.end(), CompNodeInfo<SPACE_DIM>());
            queue.pop_back();
        }
        
    }
    
};


#endif //_DECIMATOR_HPP_
