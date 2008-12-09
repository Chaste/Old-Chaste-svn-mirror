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
     }
     
    void TestVertexElementAreaAndPerimeterOnCircle()
    {
        std::vector<Node<2>*> nodes;
        unsigned N = 1000;   //vertices
        for(unsigned i=0; i<N; i++)
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
     }
     
};
#endif /*TESTVERTEXELEMENT_HPP_*/
