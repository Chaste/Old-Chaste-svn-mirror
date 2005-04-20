#ifndef _TESTNODE_HPP_
#define _TESTNODE_HPP_

#include "Point.hpp"
#include "Node.hpp"

class TestNode : public CxxTest::TestSuite 
{
	public:
	/**
	 * Test that get and set methods work as expected.
	 * Also check constructors.
	 * We only test 1 dimensional nodes. Nothing much changes in higher
	 * dimensions.
	 */
	void testNode(void)
	{
		Point<1> point1(1.0);
		Point<1> point2(2.0);
		
		Node<1> node1(0, point1);
		TS_ASSERT_EQUALS(node1.GetIndex(), 0);
		TS_ASSERT_DELTA(1.0, node1.GetPoint()[0], 1e-12);
		
		node1.SetIndex(1);
		TS_ASSERT_EQUALS(node1.GetIndex(), 1);
		
		node1.SetPoint(point2);
		TS_ASSERT_DELTA(2.0, node1.GetPoint()[0], 1e-12);
		
		node1.SetAsBoundaryNode();
		TS_ASSERT(node1.IsBoundaryNode());
		
		node1.SetAsBoundaryNode(false);
		TS_ASSERT(!node1.IsBoundaryNode());
	}
};

#endif //_TESTNODE_HPP_
