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
#include "RandomNumberGenerator.hpp"

#include <vector>
#include <cmath>

//double ComputePolarAngle(double x, double y)
//{
//    if (x==0)
//    {
//        if (y>0)
//        {
//            return M_PI/2.0;
//        }
//        else if (y<0)
//        {
//            return -M_PI/2.0;
//        }
//        else
//        {
//            EXCEPTION("Tried to compute polar angle of (0,0)");
//        }
//    }
//
//    double angle = atan(y/x);
//
//    if (y >= 0 && x < 0 )
//    {
//        angle += M_PI;
//    }
//    else if (y < 0 && x < 0 )
//    {
//        angle -= M_PI;
//    }
//    return angle;
//};

/**
 * When creating an element within a mesh one needs to specify its global index.
 * If the element is not used within a mesh the following constant is used instead.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexElement : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{
private:
     double mVertexElementArea;
     double mVertexElementPerimeter;

public:
    
    VertexElement(unsigned index, std::vector<Node<SPACE_DIM>*> nodes)
        : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index, nodes)
    {
        RegisterWithNodes();
        mVertexElementArea = DOUBLE_UNSET;
        mVertexElementPerimeter = DOUBLE_UNSET;
    }

    ~VertexElement()
    {}
    
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
    
    /** Deletes node from the given index
     *  @param rIndex is an local index to which node to remove
     */
    void DeleteNode(const unsigned& rIndex)
    {
        assert(rIndex < this->mNodes.size());
        
        // Remove element from the node at this location
        this->mNodes[rIndex]->RemoveElement(this->mIndex);
        
        // Remove the node at rIndex (removes node from element)
        this->mNodes.erase( this->mNodes.begin( ) + rIndex );

        // Update perimeter and area
        CalculateVertexElementAreaAndPerimeter();
  }
    
    /** Adds a node on the edge between nodes at rIndex and rIndex + 1
     *  @param rIndex is an local index after which node is added
     */
    void DivideEdge(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
    {
        assert(rIndex < this->mNodes.size());

        // Update the pNode location
        c_vector<double, SPACE_DIM> position;
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            position[i] = 0.5*(this->mNodes[rIndex]->GetPoint()[i]+this->mNodes[rIndex+1]->GetPoint()[i]);
        }
        ChastePoint<SPACE_DIM> point(position);
        
        pNode->SetPoint(point);
        
        // adds pNode to rIndex+1 element of mNodes pushing the others up.
        this->mNodes.insert( this->mNodes.begin( ) + rIndex+1,  pNode);

        // Add element to this node
        this->mNodes[rIndex+1]->AddElement(this->mIndex);
    }
    
     /**
     *  Calculate the area and perimeter of the polygon
     *  Assumes SPACE_DIM = 2
     */
    void CalculateVertexElementAreaAndPerimeter()
    {
        assert(SPACE_DIM == 2);
        c_vector<double,SPACE_DIM> current_node,anticlockwise_node; 
        
        double temp_vertex_element_area = 0;
        double temp_vertex_element_perimeter = 0;
        unsigned number_of_nodes = this->GetNumNodes();
               
        for (unsigned i=0; i<(number_of_nodes); i++)
        {
            // Find locations of current node and anticlockwise node
            current_node = this->GetNodeLocation(i);
            anticlockwise_node = this->GetNodeLocation((i+1)%number_of_nodes);
            
            /// \todo will need to change length calculation to something like GetVectorFromAtoB (see #825)
            
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
    
    /**
     *  Calculate the seconds moments of area of the polygon
     *  Assumes SPACE_DIM = 2
     *  @return (Ixx,Iyy,Ixy).
     */
    c_vector<double, 3> CalculateMoments()
    {
        assert(SPACE_DIM == 2);
        c_vector<double, 3> moments = zero_vector<double>(3);
        unsigned node_1, node_2;
        unsigned N = this->GetNumNodes();

        for (unsigned i=0; i<N; i++)
        {
            node_1 = i;
            node_2 = (i+1)%N;

            c_vector<double, 2> pos_1 = this->mNodes[node_1]->rGetLocation();
            c_vector<double, 2> pos_2 = this->mNodes[node_2]->rGetLocation();

            // Ixx
            moments(0) += (pos_2(0)-pos_1(0))*(  pos_1(1)*pos_1(1)*pos_1(1)
                                               + pos_1(1)*pos_1(1)*pos_2(1)
                                               + pos_1(1)*pos_2(1)*pos_2(1)
                                               + pos_2(1)*pos_2(1)*pos_2(1));

            // Iyy
            moments(1) += (pos_2(1)-pos_1(1))*(  pos_1(0)*pos_1(0)*pos_1(0)
                                               + pos_1(0)*pos_1(0)*pos_2(0)
                                               + pos_1(0)*pos_2(0)*pos_2(0)
                                               + pos_2(0)*pos_2(0)*pos_2(0));

            // Ixy
            moments(2) +=   pos_1(0)*pos_1(0)*pos_2(1)*(pos_1(1)*2 + pos_2(1))
                          - pos_2(0)*pos_2(0)*pos_1(1)*(pos_1(1) + pos_2(1)*2)
                          + 2*pos_1(0)*pos_2(0)*(pos_2(1)*pos_2(1) - pos_1(1)*pos_1(1));
        }

        moments(0) /= -12;
        moments(1) /= 12;
        moments(2) /= 24;

        return moments;
    }
    
    /**
     *  Calculate the centroid of the polygon
     *  Assumes SPACE_DIM = 2
     *  @return (centroid_x,centroid_y).
     */
    c_vector<double, SPACE_DIM> CalculateCentroid()
    {
        assert(SPACE_DIM == 2);
        c_vector<double, SPACE_DIM> centroid=zero_vector<double>(SPACE_DIM);
        c_vector<double, SPACE_DIM> current_node,anticlockwise_node; 
        
        double temp_centroid_x = 0;
        double temp_centroid_y = 0;
         
        unsigned number_of_nodes = this->GetNumNodes();
               
        for (unsigned i=0; i<(number_of_nodes); i++)
        {
            // Find locations of current node and anticlockwise node
            current_node = this->GetNodeLocation(i);
            anticlockwise_node = this->GetNodeLocation((i+1)%number_of_nodes);
            
            /// \todo will need to change length calculation to something like GetVectorFromAtoB (see #825)
            
            temp_centroid_x += (current_node[0]+anticlockwise_node[0])*(current_node[0]*anticlockwise_node[1]-current_node[1]*anticlockwise_node[0]);
            temp_centroid_y += (current_node[1]+anticlockwise_node[1])*(current_node[0]*anticlockwise_node[1]-current_node[1]*anticlockwise_node[0]);               
            
        }
        
        //double vertex_area = mVertexElementArea;
        double vertex_area = this->VertexElement<2u,2u>::GetVertexElementArea();
        double centroid_coefficient = 1.0/6.0/vertex_area;
        
        centroid(0) = centroid_coefficient*temp_centroid_x;
        centroid(1) = centroid_coefficient*temp_centroid_y;
        
        return centroid;
        
    }
    
    /**
     *  Calculate the vector of the shortest axis - this is the eigenvector 
     *  associated with the largest eigenvalue of the inertial tensor.
     *  If the polygon is regular then the eigenvalues are the same so we
     *  return a random unit vector. 
     *  Assumes SPACE_DIM = 2
     *  @return (short_axis_x, short_axis_y).
     */
    c_vector<double, SPACE_DIM> CalculateShortAxis()
    {
        assert(SPACE_DIM == 2);
        c_vector<double, SPACE_DIM> short_axis = zero_vector<double>(SPACE_DIM);
        
        c_vector<double, 3> moments = this->VertexElement<2u,2u>::CalculateMoments();
        
        double largest_eigenvalue, discriminant;            
        
        discriminant = sqrt((moments(0) - moments(1))*(moments(0) - moments(1)) + 4.0*moments(2)*moments(2));
        
        // This is always the largest eigenvalue as both eigenvalues are real as it is a 
        // symmetric matrix
        largest_eigenvalue = ((moments(0) + moments(1)) + discriminant)*0.5;       
                       
        if (fabs(discriminant)< 1e-10)
        {
            // Returning a random unit vector.
            short_axis(0) = RandomNumberGenerator::Instance()->ranf();
            short_axis(1) = sqrt(1.0-short_axis(0)*short_axis(0));
        }
        
        else
        {                      
            if (moments(2) == 0.0)
            {
                 short_axis(0) = 0.0;
                 short_axis(1) = 1.0;                       
            }
                         
            else
            {
                short_axis(0) = 1.0;
                short_axis(1) = (moments(0) - largest_eigenvalue)/moments(2);
                
                double length_short_axis = norm_2(short_axis);
                
                short_axis /= length_short_axis;
            }     
        }                                           
            
        return short_axis;
        
    }
    
    
};

#endif /*VERTEXELEMENT_HPP_*/
