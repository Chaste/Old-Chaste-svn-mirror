#ifndef _BOUNDARYELEMENT_HPP_
#define _BOUNDARYELEMENT_HPP_

#include "AbstractElement.cpp"

template <int ELEMENT_DIM, int SPACE_DIM>
class BoundaryElement : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{
public:
    BoundaryElement(unsigned index, std::vector<Node<SPACE_DIM>*> nodes, int orderOfBasisFunctions=1);
     /**
     * Create a new boundary element from a Node
     * The element has ELEMENT_DIM=0 and
     * SPACE_DIM identical to that of the node from which it is constructed
     * 
     */
    BoundaryElement(unsigned index, Node<SPACE_DIM> *node);
    void RegisterWithNodes();
    void MarkAsDeleted();
};



#endif //_BOUNDARYELEMENT_HPP_

