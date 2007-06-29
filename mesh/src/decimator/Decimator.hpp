#ifndef _DECIMATOR_HPP_
#define _DECIMATOR_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "NodeInfo.hpp"
#include "Exception.hpp"
#include <vector>
#include <ext/algorithm>

/*** \todo
* makeheap efficiency
* inner product interpolation
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
    
    // Save typing!
    typedef typename Node<SPACE_DIM>::ContainingElementIterator ContainingElementIterator;
    
    void Initialise(NodeInfo<SPACE_DIM> *pNodeInfo)
    {
        //Pre-requesite mQueue is such that mQueue[i]->GetIndex() == i
        Node<SPACE_DIM> *p_node=pNodeInfo->mpNode;
        for (ContainingElementIterator it = p_node->ContainingElementsBegin();
             it != p_node->ContainingElementsEnd();
             ++it)
        {
            Element <SPACE_DIM, SPACE_DIM> *p_element = mpMesh->GetElement(*it);
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
        for (ContainingElementIterator it = p_node->ContainingElementsBegin();
             it != p_node->ContainingElementsEnd();
             ++it)
        {
            Element<SPACE_DIM,SPACE_DIM> *p_element = mpMesh->GetElement(*it);
            measure += p_element->GetJacobianDeterminant();
        }
        if (before)
        {
            mMeasureBefore=measure;
        }
        else
        {
            mMeasureAfter=measure;
        }
    }
    
    double CalculateNeighbourhoodVolume(Node<SPACE_DIM> *pNode)
    {
        double measure=0.0;
        for (ContainingElementIterator it = pNode->ContainingElementsBegin();
             it != pNode->ContainingElementsEnd();
             ++it)
        {
            Element<SPACE_DIM,SPACE_DIM> *p_element = mpMesh->GetElement(*it);
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
 			//Put it back       
 			if (pNodeInfo->mpNode->IsBoundaryNode())
 			{
 				mpMesh->SetNode(pNodeInfo->GetIndex(), point);
 			}
        	//Only really necessary if their are boundary nodes involved
        	//since they (presently) have memory as to their use
        }
        //Put it back
        if (!pNodeInfo->mpNode->IsBoundaryNode())
 		{
 			mpMesh->SetNode(pNodeInfo->GetIndex(), point);
 		}
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
            //UpdateHeap(p_neighbour_info);
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
        OutputFileHandler handler(filePathName+"/vis_results/",true);
        
        mNodeFile=handler.OpenOutputFile("results.viznodes");
        mElementFile=handler.OpenOutputFile("results.vizelements");
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
    
    bool StoppingCondition()
    {
    	return ( (mQueue[0]->mScore >= mThresholdScore) || ExtraStoppingCondition() );
    }
    virtual bool ExtraStoppingCondition()
    {
     	return false;
    }
    
    void MakeHeap()
    {
        for (unsigned i=(mQueue.size())/2; i>0; i--)
        {
            Heapify(i);
        }
        Heapify(0);
    }
    
    /**Make a local alteration to the heap after the score of a node has changed 
     * If the score has increased then it should bubble down the minheap - run heapify
     * If the score has decrease then it should bubble up the minheap.
     * This will potentially break the heap invariant at the top of the heap.
     * 
     * @param The NodeInfo information of the node whose score has changed.
     */ 
    void UpdateHeap(NodeInfo<SPACE_DIM> *pNodeInfo)
    {
        unsigned index=pNodeInfo->mPositionInVector;
        //Bubble low values upwards
        index=HeapifyUpwards(index);
        //Bubble high values downwards
        Heapify(index);
    }
    
    /**Preserve heap invariant locally by bubbling low values up the tree.
     * Note that when a low value arrives at its new destination it does so 
     * because it was lower than its parent.  It shouldn't break the invariant
     * with the other child, but in the case of equality it may.  Therefore it's
     * a good idea to heapify.
     * 
     * @param index of the node to start this local operation
     * @return index of the last node affected
     */
    unsigned HeapifyUpwards(unsigned index)
    {
        return index;//todo
    }
     
       
    /** Core function of the minheap
     * Preserve heap invariant locally by bubbling high values down the tree.
     * @param index of the node to start this local operation
     */
    void Heapify(unsigned index)
    {
        //Set child to be the left child
        unsigned size=mQueue.size();
        unsigned child=(2*index+1)<size?(2*index+1):0;
        
        while (child)
        {
            unsigned right_child=child+1<size?child+1:0;
            
            //Set child to be the minimum of the two children
            if (right_child != 0)
            {
                if (!(*mQueue[right_child] > *mQueue[child]))
                {
                    child=right_child;
                }
            }
            
            //Check the heap invariant and fix as necessary
            if (*mQueue[index] > *mQueue[child])
            {
                //swap
                NodeInfo<SPACE_DIM> *temp=mQueue[index];
                mQueue[index]=mQueue[child];
                mQueue[index]->mPositionInVector=index;
                mQueue[child]=temp;
                mQueue[child]->mPositionInVector=child;
            }
            else
            {
                break;
            }
            index=child;
            child=(2*index+1)<size?(2*index+1):0;
        }
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
    
    void Decimate(bool interrogate=false)
    {
        while (!StoppingCondition())
        {
            ActivateOnce();
            if (interrogate)
            {
            	Interrogate();
            }
        }
    }
    void DecimateAnimate(std::string filePathName, unsigned step=1)
    {
        unsigned time=0;
        assert(SPACE_DIM<=2);
        OpenAnimationFiles(filePathName);
        
        WriteVisualiseFiles(time);
        while (!StoppingCondition() )
        {
            ActivateOnce();
            time++;
            if ( time%step == 0)
            {
                WriteVisualiseFiles(time);
            }
        }
        
        //Make sure that the final step is always shown
        if (step != 1 && time%step != 0)
        {
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
            NodeInfo<SPACE_DIM> *p_node_info=new NodeInfo<SPACE_DIM>(mpMesh->GetNode(i), i);
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
           
        //My make heap
        MakeHeap();
        /*
        //Check things are in the right places
        for (unsigned i=0; i<mQueue.size();i++)
        {
            assert(mQueue[i]->mPositionInVector == i);
        }

        for (unsigned i=0; i<mQueue.size();i++)
        {
            std::cout<<"i "<<i<<"("<<mQueue[i]->mScore<<")\t";
            if (2*(i+1)-1 < mQueue.size())
            {
                std::cout<<"l "<<2*(i+1)-1<<"("<<mQueue[2*(i+1)-1]->mScore<<")\t";
                assert(mQueue[i]->mScore <= mQueue[2*(i+1)-1]->mScore); //Left child
            }
            if (2*(i+1) < mQueue.size())
            {
                std::cout<<"r "<<2*(i+1)<<"("<<mQueue[2*(i+1)]->mScore<<")\t";
                assert(mQueue[i]->mScore <= mQueue[2*(i+1)]->mScore); //Right child
            }
            std::cout<<"\n";
        }*/
        
    }
    virtual ~Decimator()
    {
        for (unsigned i=0; i<mQueue.size(); i++)
        {
            delete mQueue[i];
        }
    }
    
#define COVERAGE_IGNORE
//Only used in debugging
    void Interrogate()
    {
        std::cout<<"Number of nodes in queue is now "<<mQueue.size()<<"\t";
        std::cout<<"Lowest score is "<<mQueue[0]->mScore
        	<<" at "<<  mQueue[0]->GetIndex();
        if (mQueue[0]->mScore == INFINITY)
        {
            std::cout<<"\n";
        }
        else
        {
            std::cout<<" to "<<  mQueue[0]->mpBestNeighbourNode->GetIndex() <<"\n";
        }        
    }
#undef COVERAGE_IGNORE
};


#endif //_DECIMATOR_HPP_
