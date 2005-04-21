#ifndef _TESTELEMENT_HPP_
#define _TESTELEMENT_HPP_

#include "Point.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include <cxxtest/TestSuite.h>

#include <vector>

class TestElement : public CxxTest::TestSuite 
{
private:
	/**
	 * Create a node with given index that has a point at the origin and is
	 * not on the boundary.
	 * 
	 * @param index The global index of the created node.
	 * @returns A pointer to a new 3d node.
	 */
	Node<3> *CreateZeroPointNode(int index)
	{
		Point<3> point(0, 0, 0);
		Node<3> *pnode = new Node<3>(index, point, false);
		return pnode;
	}

public:
	void TestConstruction()
	{
		std::vector<Node<3>*> nodes;
		for (int i=0; i<4; i++)
		{
			nodes.push_back(CreateZeroPointNode(i)); // TODO: Will need to be changed
		}
		Element<3,3> element(nodes, true);
		
		// Check nodes on the new element have the right indices
		for (int i=0; i<4; i++)
		{
			TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(i), i);
		}
		
		// Check lower order elements are created with the expected nodes.
		// Not sure if we really want to specify this, but it ensures nothing
		// has changed from the earlier code, just in case.
		for (int i=0; i < 4; i++)
        {
            for(int j=0; j < 3; j++)
            {
               TS_ASSERT_EQUALS(element.GetLowerOrderElement(i)->GetNodeGlobalIndex(j),(i+j+1) % 4);
               TS_ASSERT_DELTA(element.GetLowerOrderElement(i)->GetNodeLocation(j, 0), 0, 1e-12);
            }
        }
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetNodeGlobalIndex(0),2);
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetNodeGlobalIndex(1),3);
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(1)->GetNodeGlobalIndex(1),1);
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetNodeGlobalIndex(0),3);

		// Test AddNode method
		Node<3>* node2 = CreateZeroPointNode(10);
		element.AddNode(node2);
		TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(4), 10);
	}
	
	void testJacobian()
	{
		std::vector<Node<1>*> nodes;
		nodes.push_back(new Node<1>(0, false, 2.0));
		nodes.push_back(new Node<1>(1, false, 2.5));
		Element<1,1> element(nodes);
		const MatrixDouble *J = element.GetJacobian();
		TS_ASSERT_DELTA((*J)(0,0), 0.5, 1e-12);
	}
	

};

#endif //_TESTELEMENT_HPP_
