#ifndef _TESTNODE_HPP_
#define _TESTNODE_HPP_

#include "Point.hpp"
#include "Node.hpp"
#include <cxxtest/TestSuite.h>

class TestNode : public CxxTest::TestSuite
{
public:
    /**
     * Test that get and set methods work as expected.
     * Also check constructors.
     * We only test 1 dimensional nodes. Nothing much changes in higher
     * dimensions.
     */
    void TestNodeMethod(void)
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
        
        Node<2> node2(1, false, 1.0, 2.0);
        TS_ASSERT_EQUALS(node2.GetIndex(), 1);
        TS_ASSERT_DELTA(node2.GetPoint()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(node2.GetPoint()[1], 2.0, 1e-12);
        TS_ASSERT(!node2.IsBoundaryNode());
        
        // This shouldn't give an error, even though we specify too many
        // coordinates: the 3rd coord should be ignored.
        Node<2> node3(2, true, 1.0, 2.0, 3.0);
        
        // Test node deletion (from a mesh) methods
        TS_ASSERT(!node1.IsDeleted());
        node1.MarkAsDeleted();
        TS_ASSERT(node1.IsDeleted());
    }
    
    /**
     * Test that we can use both the new and old interfaces
     */
     
     void TestNodeNewAndOld(void)
    {
        Point<1> point1(1.0);
        
        // create a node with old interface
        Node<1> node1(0, point1);
        // check location with new interface
        TS_ASSERT_EQUALS(node1.rGetLocation()[0], point1[0]);
        
        c_vector<double, 1> new_location;
        new_location[0]=10.0;
        
        // update location with new interface
        node1.rGetModifiableLocation()=new_location;
        // check location with old interface
        TS_ASSERT_EQUALS(node1.GetPoint()[0], 10.0);
        
        // update location with old interface       
        node1.SetPoint(point1);
        // check location with new interface
        TS_ASSERT_EQUALS(node1.rGetLocation()[0], point1[0]);
    }
    
    void TestNodeNew(void)
    {
        c_vector<double, 1> point3;
        point3[0]=23.0;
        
        // create node
        Node<1> node1(0, point3);
        // read back location
        TS_ASSERT_EQUALS(node1.rGetLocation()[0], point3[0]);
        
        c_vector<double, 1> new_location;
        new_location[0]=10.0;
        
        // update location
        node1.rGetModifiableLocation()=new_location;
        // read back location
        TS_ASSERT_EQUALS(node1.rGetLocation()[0], 10.0);
    }
};

#endif //_TESTNODE_HPP_
