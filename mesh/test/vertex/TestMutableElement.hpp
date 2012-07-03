/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTPOTTSELEMENT_HPP_
#define TESTPOTTSELEMENT_HPP_

#include <cxxtest/TestSuite.h>

#include "MutableElement.hpp"

class TestMutableElement : public CxxTest::TestSuite
{
public:

    void Test1dMutableElement()
    {
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, true, 0.0));
        nodes.push_back(new Node<1>(1, false, 1.0));

        MutableElement<1> element(0, nodes);

        // Test RegisterWithNodes()
        element.RegisterWithNodes();
        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 1u);
        }

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), 1u);

        // Test IsElementOnBoundary()
        TS_ASSERT_EQUALS(element.IsElementOnBoundary(),true);

        // Test UpdateNode()
        Node<1>* p_node_2 = new Node<1>(2, false, 2.0);
        element.UpdateNode(0, p_node_2);

        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 2.0, 1e-12);

        // Test ResetIndex()
        TS_ASSERT_EQUALS(element.GetIndex(), 0u);
        element.ResetIndex(5);
        TS_ASSERT_EQUALS(element.GetIndex(), 5u);

        // Test DeleteNode() and AddNode()
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        element.DeleteNode(1);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 1u);

        Node<1>* p_node_3 = new Node<1>(3, false, 3.0);
        element.AddNode(p_node_3);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[0], 3.0, 1e-12);

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(2), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(3), 1u);

        // Test IsElementOnBoundary()
        TS_ASSERT_EQUALS(element.IsElementOnBoundary(),false);

        // Test MarkAsDeleted()
        element.MarkAsDeleted();

        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 0u);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        delete p_node_2;
        delete p_node_3;
    }

    void Test2dMutableElement()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.0, 1.0));

        MutableElement<2> element(0, nodes);

        // Test RegisterWithNodes()
        element.RegisterWithNodes();
        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 1u);
        }

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), 1u);

        // Test IsElementOnBoundary()
        TS_ASSERT_EQUALS(element.IsElementOnBoundary(),true);

        // Test UpdateNode()
        Node<2>* p_node_2 = new Node<2>(2, false, 1.0, 0.0);
        element.UpdateNode(0, p_node_2);

        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 0.0, 1e-12);

        // Test ResetIndex()
        TS_ASSERT_EQUALS(element.GetIndex(), 0u);
        element.ResetIndex(5);
        TS_ASSERT_EQUALS(element.GetIndex(), 5u);

        // Test DeleteNode() and AddNode()
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        element.DeleteNode(1);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 1u);

        Node<2>* p_node_3 = new Node<2>(3, false, 1.0, 1.0);
        element.AddNode(p_node_3);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[1], 1.0, 1e-12);

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(2), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(3), 1u);

        // Test MarkAsDeleted()
        element.MarkAsDeleted();

        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 0u);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        delete p_node_2;
        delete p_node_3;
    }

    void Test3dMutableElement()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 0.0, 1.0, 0.0));

        MutableElement<3> element(0, nodes);

        // Test RegisterWithNodes()
        element.RegisterWithNodes();
        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
           TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 1u);
        }

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), 1u);

        // Test IsElementOnBoundary()
        TS_ASSERT_EQUALS(element.IsElementOnBoundary(),true);

        // Test UpdateNode()
        Node<3>* p_node_2 = new Node<3>(2, false, 1.0, 0.0, 0.0);
        element.UpdateNode(0, p_node_2);

        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[2], 0.0, 1e-12);

        // Test ResetIndex()
        TS_ASSERT_EQUALS(element.GetIndex(), 0u);
        element.ResetIndex(5);
        TS_ASSERT_EQUALS(element.GetIndex(), 5u);

        // Test DeleteNode() and AddNode()
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        element.DeleteNode(1);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 1u);

        Node<3>* p_node_3 = new Node<3>(3, false, 1.0, 1.0, 0.0);
        element.AddNode(p_node_3);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[2], 0.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[1], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[2], 0.0, 1e-12);

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(2), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(3), 1u);

        // Test MarkAsDeleted()
        element.MarkAsDeleted();

        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
           TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 0u);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
           delete nodes[i];
        }
        delete p_node_2;
        delete p_node_3;
    }
};

#endif /*TESTPOTTSELEMENT_HPP_*/