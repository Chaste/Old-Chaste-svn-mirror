#include "NodeBox.hpp"

template<unsigned DIM>
NodeBox<DIM>::NodeBox(c_vector<double, 2*DIM> minAndMaxValues)
{
    mMinAndMaxValues = minAndMaxValues;
}

template<unsigned DIM>
c_vector<double, 2*DIM>& NodeBox<DIM>::rGetMinAndMaxValues()
{
    return mMinAndMaxValues;
}

template<unsigned DIM>
void NodeBox<DIM>::AddNode(Node<DIM>* p_node)
{
    mNodesContained.insert(p_node);
}
   
   
template<unsigned DIM>
void NodeBox<DIM>::RemoveNode(Node<DIM>* p_node)
{
    mNodesContained.erase(p_node);
}

template<unsigned DIM>
std::set< Node<DIM>* >& NodeBox<DIM>::rGetNodesContained()
{
    return mNodesContained;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class NodeBox<1>;
template class NodeBox<2>;
template class NodeBox<3>;
