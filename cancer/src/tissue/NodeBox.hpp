/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
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
