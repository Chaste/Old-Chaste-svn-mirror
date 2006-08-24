#ifndef _ELEMENT_CPP_
#define _ELEMENT_CPP_
#include "Element.hpp"

template<int ELEMENT_DIM,int SPACE_DIM>
Element<ELEMENT_DIM,SPACE_DIM>::Element(unsigned index,
        std::vector<Node<SPACE_DIM>*> nodes,
        int orderOfBasisFunctions)
        : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index,nodes,orderOfBasisFunctions)
{
    RegisterWithNodes();
 }
    
/***
 * Copy constructor which allows a new index to be specified
 */
template<int ELEMENT_DIM,int SPACE_DIM>
Element<ELEMENT_DIM,SPACE_DIM>::Element(const Element &element, const unsigned index)
{
    this->mIndex=index;
    CommonConstructor(element);
    RegisterWithNodes();
}

template<int ELEMENT_DIM,int SPACE_DIM>
void Element<ELEMENT_DIM,SPACE_DIM>::RegisterWithNodes()
{
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->AddElement(this->mIndex);
    }
}

template<int ELEMENT_DIM,int SPACE_DIM>
void Element<ELEMENT_DIM,SPACE_DIM>::MarkAsDeleted()
{
    this->mIsDeleted = true;
    this->mJacobianDeterminant = 0.0;
    // Update nodes in this element so they know they are not contained by us
    for (int i=0; i<this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveElement(this->mIndex);
    }
}

#endif //_ELEMENT_CPP_
