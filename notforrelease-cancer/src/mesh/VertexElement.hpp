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
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexElement : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{
private:

     /** Area of the element. */
     double mVertexElementArea;
     
     /** Perimeter of the element. */
     double mVertexElementPerimeter;

public:

    /**
     * Constructor.
     * 
     * @param index global index of the element
     * @param nodes vector of Nodes associated with the element
     */
    VertexElement(unsigned index, std::vector<Node<SPACE_DIM>*> nodes);

    /**
     * Destructor.
     */
    ~VertexElement();

    /**
     * Overridden RegisterWithNodes() method.
     * 
     * Informs all nodes forming this element that they are in this element.
     */
    void RegisterWithNodes();
    
    /** 
     * \todo should this return a reference? (see #861)
     * 
     * @return mNodes 
     */
    std::vector<Node<SPACE_DIM>*> GetNodes();

    /**
     * Overridden MarkAsDeleted() method.
     * 
     * Mark an element as having been removed from the mesh.
     * Also notify nodes in the element that it has been removed.
     */
    void MarkAsDeleted();

    /** 
     * Update node at the given index.
     * 
     * @param rIndex is an local index to which node to change
     * @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);

    /** 
     * Delete a node with given local index.
     * 
     * @param rIndex is the local index of the node to remove
     */
    void DeleteNode(const unsigned& rIndex);

    /** 
     * Add a node on the edge between nodes at rIndex and rIndex+1.
     * 
     * @param rIndex the local index of the node after which the new node is added
     * @param pNode a pointer to the new node
     */
    void DivideEdge(const unsigned& rIndex, Node<SPACE_DIM>* pNode);

    /**
     * Calculate the area and perimeter of the (polygonal) element, 
     * and store as the member variables mVertexElementArea and 
     * mVertexElementPerimeter.
     * 
     * \todo This method currently assumes SPACE_DIM = 2 (see #825)
     */
    void CalculateVertexElementAreaAndPerimeter();

    /**
     * @return mVertexElementArea.
     */
    double GetArea();

    /**
     * @return mVertexElementPerimeter.
     */        
    double GetPerimeter();

    /**
     * Compute the second moments of area of the (polygonal) element.
     * 
     * \todo This method currently assumes SPACE_DIM = 2 (see #825)
     * 
     * @return (Ixx,Iyy,Ixy).
     */
    c_vector<double, 3> CalculateMoments();

    /**
     * Compute the centroid of the (polygonal) element.
     * 
     * \todo This method currently assumes SPACE_DIM = 2 (see #825)
     * 
     * @return (centroid_x,centroid_y).
     */
    c_vector<double, SPACE_DIM> CalculateCentroid();

    /**
     * Calculate the vector of the shortest axis of the element. 
     * This is the eigenvector associated with the largest eigenvalue 
     * of the inertial tensor. If the polygon is regular then the 
     * eigenvalues are the same, so we return a random unit vector.
     *  
     * \todo This method currently assumes SPACE_DIM = 2 (see #825)
     *
     *  @return (short_axis_x, short_axis_y).
     */
    c_vector<double, SPACE_DIM> CalculateShortAxis();

};

#endif /*VERTEXELEMENT_HPP_*/
