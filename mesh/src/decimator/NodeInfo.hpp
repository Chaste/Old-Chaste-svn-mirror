#ifndef _NODEINFO_HPP_
#define _NODEINFO_HPP_


#include "Node.hpp"

//const double INFINITY=DBL_MAX;

template <unsigned SPACE_DIM>
class Decimator;
template <unsigned SPACE_DIM>
class CompNodeInfo;

template <unsigned SPACE_DIM>
class NodeInfo
{
private:
    friend class Decimator<SPACE_DIM>;
    friend class CompNodeInfo<SPACE_DIM>;
protected:
    double mScore;
    double mNeighbourhoodVolume;
    double mMeasureBefore;
    unsigned mPossibleTargetIndex;
    NodeInfo<SPACE_DIM> *mpBestNeighbourNode;
    Node<SPACE_DIM> *mpNode;
    std::set<NodeInfo<SPACE_DIM> *> mNeighbourNodes;
    typedef typename std::set<NodeInfo<SPACE_DIM> *>::const_iterator NeighbourNodeIterator;
    NeighbourNodeIterator mNeighbourNodeIterator;
    std::set<unsigned> mNeighbourElements;
    std::set<unsigned>::const_iterator mNeighbourElementIterator;
    
public:
    NodeInfo(Node<SPACE_DIM> *pNode) : mScore(INFINITY), mpBestNeighbourNode(0),mpNode(pNode)
{}
    
    unsigned GetPossibleTargetIndex()
    {
        return mPossibleTargetIndex;
    }
    double GetNeighbourhoodVolume()
    {
        return mNeighbourhoodVolume;
    }
    
    Node<SPACE_DIM> *pGetNode()
    {
        return mpNode;
    }
    
    void AddToNeighbourNodes(NodeInfo<SPACE_DIM> *p_neighbour)
    {
        if (p_neighbour == this)
        {
            return;
        }
        mNeighbourNodes.insert(p_neighbour);
        mNeighbourNodeIterator = mNeighbourNodes.begin();
    }
    void RemoveFromNeighbourNodes(NodeInfo<SPACE_DIM> *p_neighbour)
    {
        mNeighbourNodes.erase(p_neighbour);
        mNeighbourNodeIterator = mNeighbourNodes.begin();//Perhaps?
    }
    
    const unsigned GetNumNeighbourNodes()
    {
        return mNeighbourNodes.size();
    }
    
    
    NodeInfo<SPACE_DIM> * GetNextNeighbourNode()
    {
        NodeInfo<SPACE_DIM> *current_neighbour = *mNeighbourNodeIterator;
        mNeighbourNodeIterator++;
        
        if (mNeighbourNodeIterator == mNeighbourNodes.end())
        {
            mNeighbourNodeIterator = mNeighbourNodes.begin();
        }
        return current_neighbour;
    }
    
    /*
    void AddToNeighbourElements(unsigned neighbour)
    {
        mNeighbourElements.insert(neighbour);
        mNeighbourElementIterator = mNeighbourElements.begin();
    }
    void RemoveFromNeighbourElements(unsigned neighbour)
    {
        mNeighbourElements.erase(neighbour);
        mNeighbourElementIterator = mNeighbourElements.begin();//Perhaps?
    }
    
    const unsigned GetNumNeighbourElements()
    {
        return mNeighbourElements.size();
    }
    
    unsigned GetNextNeighbourElement()
    {
        unsigned current_neighbour = *mNeighbourElementIterator;
        mNeighbourElementIterator++;
        
        if (mNeighbourElementIterator == mNeighbourElements.end())
        {
            mNeighbourElementIterator = mNeighbourElements.begin();
        }
        return current_neighbour;
    }
    */
    
    unsigned GetIndex()
    {
        return mpNode->GetIndex();
    }
    
    double GetScore()
    {
        return mScore;
    }
};

template <unsigned SPACE_DIM>
class CompNodeInfo : public std::binary_function<NodeInfo<SPACE_DIM>*, NodeInfo<SPACE_DIM> *, bool>
{
public:
    bool operator () (const NodeInfo<SPACE_DIM>* p1, const NodeInfo<SPACE_DIM>* p2)
    {
        return p1->mScore > p2->mScore;
    }
};


#endif //_NODEINFO_HPP_
