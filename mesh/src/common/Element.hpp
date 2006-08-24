#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

#include "AbstractElement.cpp"

template <int ELEMENT_DIM, int SPACE_DIM>
class Element : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{
 public:
    Element(unsigned index, std::vector<Node<SPACE_DIM>*> nodes, int orderOfBasisFunctions=1);
    /***
     * Copy constructor which allows a new index to be specified
     */
    Element(const Element &element, const unsigned index);
    void RegisterWithNodes();
    void MarkAsDeleted();
 };



#endif //_BOUNDARYELEMENT_HPP_

