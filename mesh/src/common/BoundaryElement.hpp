#ifndef _BOUNDARYELEMENT_HPP_
#define _BOUNDARYELEMENT_HPP_

#include "AbstractElement.cpp"

template <int ELEMENT_DIM, int SPACE_DIM>
class BoundaryElement : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{

public:
    BoundaryElement(unsigned index,
                    std::vector<Node<SPACE_DIM>*> nodes,
                    int orderOfBasisFunctions=1): AbstractElement<ELEMENT_DIM, SPACE_DIM>(index,nodes,orderOfBasisFunctions)
    {
        RegisterWithNodes();
    }
    
    /**
     * Create a new boundary element from a Node
     * The element has ELEMENT_DIM=0 and
     * SPACE_DIM identical to that of the node from which it is constructed
     * 
     */
    BoundaryElement(unsigned index,
                    Node<SPACE_DIM> *node)
    {
        assert (ELEMENT_DIM == 0);
        
        // Store Node pointer
        this->mNodes.push_back(node);
        
        this->mpJacobian = NULL;
        this->mpInverseJacobian = NULL;
        this->mpJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
        this->mpInverseJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
        (*this->mpJacobian)(0,0) = 1.0;
        (*this->mpInverseJacobian)(0,0) = 1.0;
        this->mJacobianDeterminant = 1.0;
        
        RegisterWithNodes();
        
    }
    
    void RegisterWithNodes()
    {
        for (unsigned i=0; i<this->mNodes.size(); i++)
        {
            this->mNodes[i]->AddBoundaryElement(this->mIndex);
        }
    }
    
    void MarkAsDeleted()
    {
        this->mIsDeleted = true;
        // Update nodes in this element so they know they are not contained by us
        for (int i=0; i<this->GetNumNodes(); i++)
        {
            this->mNodes[i]->RemoveBoundaryElement(this->mIndex);
        }
    }
};



#endif //_BOUNDARYELEMENT_HPP_

