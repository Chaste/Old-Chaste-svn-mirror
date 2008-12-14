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

 void TestVertexElementDeleteNode()
    {
//        std::vector<Node<2>*> corner_nodes;
//        corner_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
//        corner_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
//        corner_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
//        corner_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
//    
//        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, corner_nodes);
//        
//        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 4u);
//        TS_ASSERT_DELTA(vertex_element.GetVertexElementArea(),1.0,1e-6);
//        TS_ASSERT_DELTA(vertex_element.GetVertexElementPerimeter(),4.0,1e-6);
//       
        std::vector<Node<2>*> nodes;
        unsigned N = 6;   //vertices
        for(unsigned i=0; i<N; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(N); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, nodes);
        
        TS_ASSERT_DELTA(vertex_element.GetVertexElementArea(),3*sqrt(3)/2.0,1e-4);
        TS_ASSERT_DELTA(vertex_element.GetVertexElementPerimeter(),6.0,1e-4);
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

        // Tests Areas updated
        TS_ASSERT_DELTA(vertex_element.GetVertexElementArea(),sqrt(3.0),1e-6);
        TS_ASSERT_DELTA(vertex_element.GetVertexElementPerimeter(),2.0+2.0*sqrt(3.0),1e-6);
                 
        for (unsigned i=0; i<nodes.size(); ++i)
        {
            delete nodes[i];
        }
     }


void TestVertexElementDivideEdge()
    {
        std::vector<Node<2>*> corner_nodes;
        corner_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        corner_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        corner_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        corner_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
    
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, corner_nodes);
        
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 4u);
        // Pass pointer to new node position of node is not important as it will be changed
        Node<2>* new_node= new Node<2>(4, false, 0.0, 0.0);
        vertex_element.DivideEdge(0, new_node); // Divide edge between nodes (0,0) and (1,0)
        
        // Test edge is divided
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 5u);
        
        TS_ASSERT_DELTA(vertex_element.GetVertexElementArea(),1.0,1e-6);
        TS_ASSERT_DELTA(vertex_element.GetVertexElementPerimeter(),4.0,1e-6);
                 
        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 0.0, 1e-9);
        
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 0.0, 1e-9);
            
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], 0.0, 1e-9);     
                 
        for (unsigned i=0; i<corner_nodes.size(); ++i)
        {
            delete corner_nodes[i];
        }
        delete new_node;
     }
     
    void TestVertexElementAreaAndPerimeter()
    {
        std::vector<Node<2>*> corner_nodes;
        corner_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        corner_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        corner_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        corner_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
    
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, corner_nodes);
        
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 4u);
        //  Check nodes have correct indices
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNodeGlobalIndex(i), i);
        }
        
        TS_ASSERT_DELTA(vertex_element.GetVertexElementArea(),1.0,1e-6);
        TS_ASSERT_DELTA(vertex_element.GetVertexElementPerimeter(),4.0,1e-6);
        
        for (unsigned i=0; i<corner_nodes.size(); ++i)
        {
            delete corner_nodes[i];
        }
     }
     
    void TestVertexElementAreaAndPerimeterOnCircle()
    {
        std::vector<Node<2>*> nodes;
        unsigned N = 1000;   //vertices
        for (unsigned i=0; i<N; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(N); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, nodes);
        
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), N);
        //  Check nodes have correct indices
        for (unsigned i=0; i<N; i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNodeGlobalIndex(i), i);
        }
        
        TS_ASSERT_DELTA(vertex_element.GetVertexElementArea(),M_PI,1e-4);
        TS_ASSERT_DELTA(vertex_element.GetVertexElementPerimeter(),2.0*M_PI,1e-4);
        
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
        std::vector<Node<2>*> nodes;
        unsigned N = 6;   //vertices
        for (unsigned i=0; i<N; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(N); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }
        VertexElement<2,2> hexagon(INDEX_IS_NOT_USED, nodes);
        
        TS_ASSERT_DELTA(hexagon.GetVertexElementArea(),3*sqrt(3)/2.0,1e-4);
        TS_ASSERT_DELTA(hexagon.GetVertexElementPerimeter(),6.0,1e-4);
        
        c_vector<double, 3> moments = hexagon.CalculateMoments();
        TS_ASSERT_DELTA(moments(0), 5*sqrt(3)/16, 1e-6);    // Ixx
        TS_ASSERT_DELTA(moments(1), 5*sqrt(3)/16, 1e-6);    // Iyy
        TS_ASSERT_DELTA(moments(2), 0.0, 1e-6);    // Ixy
        
        for (unsigned i=0; i<nodes.size(); ++i)
        {
            delete nodes[i];
        }
    }
    
    void TestCalculateCentroid() throw(Exception)
    {
        std::vector<Node<2>*> nodes;
        unsigned N = 6;   //vertices
        for(unsigned i=0; i<N; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(N); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        VertexElement<2,2> hexagon(INDEX_IS_NOT_USED, nodes);
       
        c_vector<double, 2> centroid = hexagon.CalculateCentroid();
        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(centroid(1), 0.0, 1e-6);
        
        for(unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
     
         void TestCalculateShortAxis() throw(Exception)
    {
        std::vector<Node<2>*> nodes1;
       
        // This is a rectangle, centre (0,0), width 2, length 2, parallel to x axis                       
        nodes1.push_back(new Node<2>(0, false,  2.0,  1.0));
        nodes1.push_back(new Node<2>(1, false, -2.0,  1.0));
        nodes1.push_back(new Node<2>(2, false, -2.0, -1.0));
        nodes1.push_back(new Node<2>(3, false,  2.0, -1.0));
        
        VertexElement<2,2> rectangle1(INDEX_IS_NOT_USED, nodes1);
       
        TS_ASSERT_DELTA(rectangle1.GetVertexElementArea(),8.0,1e-4);
        TS_ASSERT_DELTA(rectangle1.GetVertexElementPerimeter(),12.0,1e-4);
        
        c_vector<double, 2> centroid = rectangle1.CalculateCentroid();
        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(centroid(1), 0.0, 1e-6); 
        
        c_vector<double, 2> short_axis = rectangle1.CalculateShortAxis();
        TS_ASSERT_DELTA(short_axis(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(short_axis(1), 1.0, 1e-6);
       
        for (unsigned i=0; i<nodes1.size(); ++i)
        {
            delete nodes1[i];
        }
       
       
        std::vector<Node<2>*> nodes2;                      
        // This is a rectangle, centre (0,0), width 1, length sqrt(3), rotated by 30 degrees anticlockwise                      
        nodes2.push_back(new Node<2>(0, false,  1.0, 0.0));
        nodes2.push_back(new Node<2>(1, false,  0.5, sqrt(3.0)/2.0));
        nodes2.push_back(new Node<2>(2, false, -1.0, 0.0));
        nodes2.push_back(new Node<2>(3, false, -0.5, -sqrt(3.0)/2.0));
        
        VertexElement<2,2> rectangle2(INDEX_IS_NOT_USED, nodes2);
       
        TS_ASSERT_DELTA(rectangle2.GetVertexElementArea(),sqrt(3.0),1e-4);
        TS_ASSERT_DELTA(rectangle2.GetVertexElementPerimeter(),2.0*sqrt(3.0)+2.0,1e-4);
        
        centroid = rectangle2.CalculateCentroid();
        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(centroid(1), 0.0, 1e-6); 
        
        short_axis = rectangle2.CalculateShortAxis();
        TS_ASSERT_DELTA(short_axis(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(short_axis(1), -sqrt(3.0)*0.5, 1e-6);
        
        for (unsigned i=0; i<nodes2.size(); ++i)
        {
            delete nodes2[i];
        }
        
        // Test on a regular polygon (generates a random vector)
        std::vector<Node<2>*> hexagon_nodes;
        unsigned N = 6;   //vertices
        for(unsigned i=0; i<N; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(N); 
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }
        VertexElement<2,2> hexagon(INDEX_IS_NOT_USED, hexagon_nodes);
        
        short_axis = hexagon.CalculateShortAxis();
        
        TS_ASSERT_DELTA(short_axis(0)*short_axis(0)+short_axis(1)*short_axis(1), 1.0, 1e-6);
        
        for (unsigned i=0; i<hexagon_nodes.size(); ++i)
        {
            delete hexagon_nodes[i];
        }
    }
};
#endif /*TESTVERTEXELEMENT_HPP_*/
