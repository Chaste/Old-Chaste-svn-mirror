/*

Copyright (C) University of Oxford, 2008

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


#ifndef VERTEXELEMENT_HPP_
#define VERTEXELEMENT_HPP_

#include "Node.hpp"
#include "ChastePoint.hpp"
#include "UblasCustomFunctions.hpp"
#include "AbstractElement.hpp"
#include "Exception.hpp"

#include <vector>
#include <cmath>

// When creating an element within a mesh one needs to specify its global index
// If the element is not used within a mesh the following
// constant is used instead.

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexElement : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{
//private:
//    unsigned mIndex;
//    std::vector<Node<SPACE_DIM>*> mNodes;
//    
//    
private:
     double mVertexElementArea;
     double mVertexElementPerimeter;

public:
    
//    ///Main constructor
//    Element(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);
    VertexElement(unsigned index, std::vector<Node<SPACE_DIM>*> nodes)
        : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index, nodes)
    {
        RegisterWithNodes();
        mVertexElementArea = DOUBLE_UNSET;
        mVertexElementPerimeter = DOUBLE_UNSET;
        
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
        for (unsigned i=0; i<this->GetNumNodes(); i++)
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
    
    void CalculateVertexElementAreaAndPerimeter()
    {
        c_vector<double,SPACE_DIM> current_node,anticlockwise_node; 
        
        double temp_vertex_element_area = 0;
        double temp_vertex_element_perimeter = 0;
        unsigned number_of_nodes = this->GetNumNodes();
               
        for (unsigned i=0; i<(number_of_nodes); i++)
        {
            // Find locations of current node and anticlockwise node
            current_node = this->GetNodeLocation(i);
            anticlockwise_node = this->GetNodeLocation((i+1)%number_of_nodes);
            
            // TODO: will need to change length calculation to something like GetVectorFromAtoB
            
            temp_vertex_element_area += 0.5*(current_node[0]*anticlockwise_node[1]
                -anticlockwise_node[0]*current_node[1]);
                
            temp_vertex_element_perimeter += norm_2(current_node-anticlockwise_node);
        }
        
        mVertexElementArea = temp_vertex_element_area;
        mVertexElementPerimeter = temp_vertex_element_perimeter;
    }
    
    double GetVertexElementArea()
    {
        if (mVertexElementArea == DOUBLE_UNSET)
        {
            this->CalculateVertexElementAreaAndPerimeter();
        }
        
        return mVertexElementArea;
    }
        
    double GetVertexElementPerimeter()
    {
        if (mVertexElementPerimeter == DOUBLE_UNSET)
        {
            this->CalculateVertexElementAreaAndPerimeter();
        }
        
        return mVertexElementPerimeter;
    }
    
};

#endif /*VERTEXELEMENT_HPP_*/
