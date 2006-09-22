#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

#include "AbstractElement.cpp"

template <int ELEMENT_DIM, int SPACE_DIM>
class Element : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{

public:
    Element(unsigned index,
            std::vector<Node<SPACE_DIM>*> nodes,
            int orderOfBasisFunctions=1): AbstractElement<ELEMENT_DIM, SPACE_DIM>(index,nodes,orderOfBasisFunctions)
    {
        RegisterWithNodes();
    }
    
    /***
     * Copy constructor which allows a new index to be specified
     */
    Element(const Element &element, const unsigned index)
    {
        this->mIndex=index;
        CommonConstructor(element);
        RegisterWithNodes();
    }
    
    void RegisterWithNodes()
    {
        for (unsigned i=0; i<this->mNodes.size(); i++)
        {
            this->mNodes[i]->AddElement(this->mIndex);
        }
    }
    
    void MarkAsDeleted()
    {
        this->mIsDeleted = true;
        this->mJacobianDeterminant = 0.0;
        // Update nodes in this element so they know they are not contained by us
        for (int i=0; i<this->GetNumNodes(); i++)
        {
            this->mNodes[i]->RemoveElement(this->mIndex);
        }
    }
    
    /** Update node at the given index
     *  @param rIndex is an local index to which node to change
     *  @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
    {
        assert(rIndex < this->mNodes.size());
        
        // Remove it from the node at this location
        this->mNodes[rIndex]->RemoveElement(this->mIndex);
        
        // Update the node at this location
        this->mNodes[rIndex] = pNode;
        
        // Add element to this node
        this->mNodes[rIndex]->AddElement(this->mIndex);
    }
    
    void ResetIndex(int index)
    {
        for (int i=0; i<this->GetNumNodes(); i++)
        {
            this->mNodes[i]->RemoveElement(this->mIndex);
        }
        this->mIndex=index;
        RegisterWithNodes();
    }
    
    double CalculateCircumsphereVolume()
    {
    
    
        return 0;
    }
    
};



#endif //_BOUNDARYELEMENT_HPP_

