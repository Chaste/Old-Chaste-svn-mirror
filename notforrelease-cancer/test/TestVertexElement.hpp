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

#ifndef TESTVERTEXELEMENT_HPP_
#define TESTVERTEXELEMENT_HPP_

#include <cxxtest/TestSuite.h>
#include "VertexElement.hpp"


class TestVertexElement : public CxxTest::TestSuite
{
public:

    void TestVertexElementDeleteAndAddNode()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }

        // Create element        
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, nodes);
        
        TS_ASSERT_DELTA(vertex_element.GetArea(), 3*sqrt(3)/2.0, 1e-4);
        TS_ASSERT_DELTA(vertex_element.GetPerimeter(), 6.0, 1e-4);
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 6u);
             
        vertex_element.DeleteNode(3); // Removes (-1,0) node
        vertex_element.DeleteNode(0); // Removes (1,0) node
        
        // Test node is removed
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 4u);
               
        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);
 
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 0.5, 1e-9);        
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        // Tests areas updated
        TS_ASSERT_DELTA(vertex_element.GetArea(), sqrt(3.0), 1e-6);
        TS_ASSERT_DELTA(vertex_element.GetPerimeter(), 2.0+2.0*sqrt(3.0), 1e-6);

        // Add new node
        Node<2>* p_new_node = new Node<2>(4, false, 0.0, 0.0);
        vertex_element.AddNode(3, p_new_node); // Add node at (0,0) between nodes 3 and 0
        
        // Test node is added
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 5u);
               
        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);
 
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 0.5, 1e-9);        
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[0], 0.0, 1e-9);        
        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[1], 0.0, 1e-9);

        // Test element area and perimeter are updated
        TS_ASSERT_DELTA(vertex_element.GetArea(), sqrt(3.0)*3.0/4.0, 1e-6);
        TS_ASSERT_DELTA(vertex_element.GetPerimeter(), 2.0+sqrt(3.0)+2.0, 1e-6);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        delete p_new_node;
    }


    void TestVertexElementDivideEdge()
    {
        // Create nodes
        std::vector<Node<2>*> corner_nodes;
        corner_nodes.push_back(new Node<2>(0, false, 1.0, 1.0));
        corner_nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        corner_nodes.push_back(new Node<2>(2, false, 2.0, 2.0));
        corner_nodes.push_back(new Node<2>(3, false, 1.0, 2.0));

        // Create element
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, corner_nodes);
        
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 4u);
        
        // Divide an edge

        // Pass pointer to new node (position of node is not important as it will be changed)
        Node<2>* p_new_node1 = new Node<2>(4, false, 0.0, 0.0);        
        vertex_element.DivideEdge(0, p_new_node1); // Divide edge between nodes (1,1) and (2,1)
        
        // Test edge is divided
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 5u);        
        TS_ASSERT_DELTA(vertex_element.GetArea(), 1.0, 1e-6);
        TS_ASSERT_DELTA(vertex_element.GetPerimeter(), 4.0, 1e-6);
        
        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 1.0, 1e-9);        
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], 1.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 1.0, 1e-9);            
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[1], 2.0, 1e-9);
        
        // Divide another edge

        Node<2>* p_new_node2 = new Node<2>(5, false, 0.0, 0.0);
        vertex_element.DivideEdge(4, p_new_node2); // Divide edge between nodes (1,2) and (1,1)
        
        // Test edge is divided
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 6u);        
        TS_ASSERT_DELTA(vertex_element.GetArea(), 1.0, 1e-6);
        TS_ASSERT_DELTA(vertex_element.GetPerimeter(), 4.0, 1e-6);
                 
        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], 1.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(5)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(5)->GetPoint()[1], 1.5, 1e-9);    

        // Tidy up
        for (unsigned i=0; i<corner_nodes.size(); i++)
        {
            delete corner_nodes[i];
        }
        delete p_new_node1;
        delete p_new_node2;
    }


    void TestVertexElementAreaAndPerimeter()
    {
        // Create nodes
        std::vector<Node<2>*> corner_nodes;
        corner_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        corner_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        corner_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        corner_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
    
        // Create element
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, corner_nodes);
        
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 4u);

        // Check nodes have correct indices
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNodeGlobalIndex(i), i);
        }
        
        // Test area and perimeter calculations
        TS_ASSERT_DELTA(vertex_element.GetArea(), 1.0, 1e-6);
        TS_ASSERT_DELTA(vertex_element.GetPerimeter(), 4.0, 1e-6);
        
        // Tidy up
        for (unsigned i=0; i<corner_nodes.size(); ++i)
        {
            delete corner_nodes[i];
        }
    }


    void TestGetAreaGradientAtNode()
    {
        // Create nodes
        std::vector<Node<2>*> corner_nodes;
        corner_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        corner_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        corner_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        corner_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, corner_nodes);
        
        // Test gradient of area evaluated at each node

        c_vector<double, 2> element_area_gradient = vertex_element.GetAreaGradientAtNode(0);        
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);
        
        element_area_gradient = vertex_element.GetAreaGradientAtNode(1);        
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);
        
        element_area_gradient = vertex_element.GetAreaGradientAtNode(2);        
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);
        
        element_area_gradient = vertex_element.GetAreaGradientAtNode(3);        
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);
        
        // Tidy up
        for (unsigned i=0; i<corner_nodes.size(); ++i)
        {
            delete corner_nodes[i];
        }
    }


    void TestVertexElementAreaAndPerimeterOnCircle()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 1000;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }

        // Create element
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, nodes);

        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), num_nodes);
        
        //  Check nodes have correct indices
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNodeGlobalIndex(i), i);
        }

        // Test area and perimeter calculations
        TS_ASSERT_DELTA(vertex_element.GetArea(), M_PI, 1e-4);
        TS_ASSERT_DELTA(vertex_element.GetPerimeter(), 2.0*M_PI, 1e-4);
        
        // Tidy up
        for (unsigned i=0; i<nodes.size(); ++i)
        {
            delete nodes[i];
        }
    }
     
//     void xTestAnticlockwisenessOfNodes() throw(Exception)
//     {
//        // Tests to check that the nodes are anticlockwise when we create element
//        std::vector<Node<2>*> corner_nodes;
//        corner_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
//        corner_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
//        corner_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
//        corner_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
//    
//        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, corner_nodes);
//        
//        std::vector<Node<2>*> corner_nodes2;
//        corner_nodes2.push_back(new Node<2>(0, false, 0.0, 0.0));
//        corner_nodes2.push_back(new Node<2>(2, false, 1.0, 1.0));   
//        corner_nodes2.push_back(new Node<2>(1, false, 1.0, 0.0));
//        corner_nodes2.push_back(new Node<2>(3, false, 0.0, 1.0));
//    
//        VertexElement<2,2> vertex_element2(INDEX_IS_NOT_USED, corner_nodes2);
//
//        for (unsigned i=0; i<corner_nodes.size(); ++i)
//        {
//            delete corner_nodes[i];
//        }
//        for (unsigned i=0; i<corner_nodes2.size(); ++i)
//        {
//            delete corner_nodes2[i];
//        }
//     }

    /// \todo this should be for a non regular polygon so Ixy ~= 0 (see #825)
    void TestCalculateMoment() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }

        // Create element
        VertexElement<2,2> hexagon(INDEX_IS_NOT_USED, nodes);
        
        TS_ASSERT_DELTA(hexagon.GetArea(),3*sqrt(3)/2.0,1e-4);
        TS_ASSERT_DELTA(hexagon.GetPerimeter(),6.0,1e-4);
        
        c_vector<double, 3> moments = hexagon.CalculateMoments();
        TS_ASSERT_DELTA(moments(0), 5*sqrt(3)/16, 1e-6);    // Ixx
        TS_ASSERT_DELTA(moments(1), 5*sqrt(3)/16, 1e-6);    // Iyy
        TS_ASSERT_DELTA(moments(2), 0.0, 1e-6);    // Ixy

        // Tidy up
        for (unsigned i=0; i<nodes.size(); ++i)
        {
            delete nodes[i];
        }
    }
    
    void TestCalculateCentroid() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        // Create element
        VertexElement<2,2> hexagon(INDEX_IS_NOT_USED, nodes);
       
        c_vector<double, 2> centroid = hexagon.CalculateCentroid();
        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(centroid(1), 0.0, 1e-6);
        
        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
     
    void TestCalculateShortAxis() throw(Exception)
    {
        // First test
        
        // Create nodes
        std::vector<Node<2>*> nodes1;
       
        // This is a rectangle, centre (0,0), width 2, length 2, parallel to x axis                       
        nodes1.push_back(new Node<2>(0, false,  2.0,  1.0));
        nodes1.push_back(new Node<2>(1, false, -2.0,  1.0));
        nodes1.push_back(new Node<2>(2, false, -2.0, -1.0));
        nodes1.push_back(new Node<2>(3, false,  2.0, -1.0));

        // Create element
        VertexElement<2,2> rectangle1(INDEX_IS_NOT_USED, nodes1);
       
        TS_ASSERT_DELTA(rectangle1.GetArea(),8.0,1e-4);
        TS_ASSERT_DELTA(rectangle1.GetPerimeter(),12.0,1e-4);

        // Test centroid calculation         
        c_vector<double, 2> centroid = rectangle1.CalculateCentroid();
        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(centroid(1), 0.0, 1e-6);

        // Test short axis calculation       
        c_vector<double, 2> short_axis = rectangle1.CalculateShortAxis();
        TS_ASSERT_DELTA(short_axis(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(short_axis(1), 1.0, 1e-6);

        // Tidy up
        for (unsigned i=0; i<nodes1.size(); ++i)
        {
            delete nodes1[i];
        }

        // Second test
        
        // Create nodes
        std::vector<Node<2>*> nodes2;
        // This is a rectangle, centre (0,0), width 1, length sqrt(3), rotated by 30 degrees anticlockwise                      
        nodes2.push_back(new Node<2>(0, false,  1.0, 0.0));
        nodes2.push_back(new Node<2>(1, false,  0.5, sqrt(3.0)/2.0));
        nodes2.push_back(new Node<2>(2, false, -1.0, 0.0));
        nodes2.push_back(new Node<2>(3, false, -0.5, -sqrt(3.0)/2.0));

        // Create element
        VertexElement<2,2> rectangle2(INDEX_IS_NOT_USED, nodes2);
       
        TS_ASSERT_DELTA(rectangle2.GetArea(),sqrt(3.0),1e-4);
        TS_ASSERT_DELTA(rectangle2.GetPerimeter(),2.0*sqrt(3.0)+2.0,1e-4);

        // Test centroid calculation        
        centroid = rectangle2.CalculateCentroid();
        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(centroid(1), 0.0, 1e-6); 

        // Test short axis calculation        
        short_axis = rectangle2.CalculateShortAxis();
        TS_ASSERT_DELTA(short_axis(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(short_axis(1), -sqrt(3.0)*0.5, 1e-6);

        // Tidy up
        for (unsigned i=0; i<nodes2.size(); ++i)
        {
            delete nodes2[i];
        }

        // Third test
               
        // Test on a regular polygon (generates a random vector)
        std::vector<Node<2>*> hexagon_nodes;
        unsigned num_nodes = 6;   // vertices
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes); 
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }

        // Create element
        VertexElement<2,2> hexagon(INDEX_IS_NOT_USED, hexagon_nodes);

        // Test short axis calculation        
        short_axis = hexagon.CalculateShortAxis();
        
        TS_ASSERT_DELTA(short_axis(0)*short_axis(0)+short_axis(1)*short_axis(1), 1.0, 1e-6);

        // Tidy up
        for (unsigned i=0; i<hexagon_nodes.size(); ++i)
        {
            delete hexagon_nodes[i];
        }
    }
    
    void TestGetNodeLocalIndex()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        // This is a square,                   
        nodes.push_back(new Node<2>(3, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(0, false, 0.0, 1.0));

        // Create element        
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, nodes);
        
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0),3u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(1),2u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(2),1u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(3),0u);
        
        vertex_element.DeleteNode(3); // Removes (1,1) node
        
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0),UINT_MAX);
        
        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
    
 
};
#endif /*TESTVERTEXELEMENT_HPP_*/
