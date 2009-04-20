#ifndef NODEBOX_HPP_
#define NODEBOX_HPP_

#include "Node.hpp"

//#include <boost/serialization/set.hpp>
//#include <boost/serialization/vector.hpp>

template<unsigned DIM>
class NodeBox
{
private:
    c_vector<double, 2*DIM> mMinAndMaxValues;
    std::set< Node<DIM>* > mNodesContained;
    
    
public:
    NodeBox(c_vector<double, 2*DIM> minAndMaxValues);
    
    c_vector<double, 2*DIM>& rGetMinAndMaxValues();
    
    void AddNode(Node<DIM>* p_node);
    void RemoveNode(Node<DIM>* p_node);
    
    std::set< Node<DIM>* >& rGetNodesContained();
};


//template<unsigned DIM>
//class NodeBoxCollection 
//{
//    /** A vector of boxes to store rough node positions */
//    std::vector< NodeBox<DIM> > mBoxes;
//    
//    c_vector<unsigned,DIM> mNumBoxesEachDirection;
//    
//    

#endif /*NODEBOX_HPP_*/
