#ifndef _TESTELEMENT_HPP_
#define _TESTELEMENT_HPP_

#include "Point.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include <cxxtest/TestSuite.h>
#include <iostream>

#include <vector>

class TestElement : public CxxTest::TestSuite 
{
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
	void dontTestConstructionForQuadraticBasisFunctions()
	{
		std::vector<const Node<3>*> cornerNodes;
//		for (int i=0; i<4; i++)
//		{
//			nodes2d.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0)); //nodes.push_back(CreateZeroPointNode(i)); // TODO: Will need to be changed
//		}
		cornerNodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
		cornerNodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
		cornerNodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
		cornerNodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
		Element<3,3> element(cornerNodes, 2,true);
		
		element.AddInternalNode(new Node<3>(4, false, 0.5, 0.0, 0.0));
		element.AddInternalNode(new Node<3>(5, false, 0.5, 0.5, 0.0));
		element.AddInternalNode(new Node<3>(6, false, 0.0, 0.5, 0.0));
		element.AddInternalNode(new Node<3>(7, false, 0.0, 0.0, 0.5));
		element.AddInternalNode(new Node<3>(8, false, 0.5, 0.0, 0.5));
		element.AddInternalNode(new Node<3>(9, false, 0.0, 0.5, 0.5));
		
		// Check nodes on the new element have the right indices
		for (int i=0; i<10; i++)
		{
			TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(i), i);
		}
		
		
		// Check lower order elements are created with the expected nodes.
		// Not sure if we really want to specify this, but it ensures nothing
		// has changed from the earlier code, just in case.
		for (int i=0; i < 10; i++)
        {
            for(int j=0; j < 6; j++)
            {
               TS_ASSERT_EQUALS(element.GetLowerOrderElement(i)->GetNodeGlobalIndex(j),(i+j+1) % 4);
               if((i==0 && j==0) || (i==2 && j==2)|| (i==3 && j==1))
               {
               		TS_ASSERT_DELTA(element.GetLowerOrderElement(i)->GetNodeLocation(j, 0), 1.0, 1e-12);
               }
               else
               {
               		TS_ASSERT_DELTA(element.GetLowerOrderElement(i)->GetNodeLocation(j, 0), 0.0, 1e-12);
               }
            }
        }
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetNodeGlobalIndex(0),2);
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetNodeGlobalIndex(1),3);
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(1)->GetNodeGlobalIndex(1),1);
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetNodeGlobalIndex(0),3);


		// Test AddNode method
		//Node<3>* node2 = CreateZeroPointNode(10);
		//element.AddNode(node2);
	
		//TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(4), 10);
	}
	
	void TestConstructionForLinearBasisFunctions()
	{
		std::vector<const Node<3>*> cornerNodes;
//		for (int i=0; i<4; i++)
//		{
//			nodes2d.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0)); //nodes.push_back(CreateZeroPointNode(i)); // TODO: Will need to be changed
//		}
		cornerNodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
		cornerNodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
		cornerNodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
		cornerNodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
		Element<3,3> element(cornerNodes, 1, true);
				
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
               if((i==0 && j==0) || (i==2 && j==2) || (i==3 && j==1))
               {
               		TS_ASSERT_DELTA(element.GetLowerOrderElement(i)->GetNodeLocation(j, 0), 1.0, 1e-12);
               }
               else
               {
               		TS_ASSERT_DELTA(element.GetLowerOrderElement(i)->GetNodeLocation(j, 0), 0.0, 1e-12);
               }
            }
        }
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetNodeGlobalIndex(0),2);
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetNodeGlobalIndex(1),3);
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(1)->GetNodeGlobalIndex(1),1);
        TS_ASSERT_EQUALS(element.GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetLowerOrderElement(0)->GetNodeGlobalIndex(0),3);


		// Test AddNode method
		//Node<3>* node2 = CreateZeroPointNode(10);
		//element.AddNode(node2);
	
		//TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(4), 10);
	}

	
	void testJacobian()
	{
		// 1d
		std::vector<const Node<1>*> nodes1d;
		nodes1d.push_back(new Node<1>(0, false, 2.0));
		nodes1d.push_back(new Node<1>(1, false, 2.5));
		Element<1,1> element1d(nodes1d);
		const MatrixDouble *J1d = element1d.GetJacobian();
		TS_ASSERT_DELTA((*J1d)(0,0), 0.5, 1e-12);
		
		double Det1d = element1d.GetJacobianDeterminant();
		TS_ASSERT_DELTA(Det1d, 0.5, 1e-12);
		const MatrixDouble *J1dinv = element1d.GetInverseJacobian();
		TS_ASSERT_DELTA((*J1dinv)(0,0), 2.0, 1e-12);
		
		// 2d easy
		
		std::vector<const Node<2>*> nodes2d;
		nodes2d.push_back(new Node<2>(0, false, 0.0, 0.0));
		nodes2d.push_back(new Node<2>(1, false, 1.0, 0.0));
		nodes2d.push_back(new Node<2>(2, false, 0.0, 1.0));
		Element<2,2> element2d(nodes2d);
		const MatrixDouble *J2d = element2d.GetJacobian();
		TS_ASSERT_DELTA((*J2d)(0,0), 1.0, 1e-12);
		TS_ASSERT_DELTA((*J2d)(0,1), 0.0, 1e-12);
		TS_ASSERT_DELTA((*J2d)(1,0), 0.0, 1e-12);
		TS_ASSERT_DELTA((*J2d)(1,1), 1.0, 1e-12);
		
		//2d general
		std::vector<const Node<2>*> nodes2d2;
		nodes2d2.push_back(new Node<2>(0, false, 1.0, -2.0));
		nodes2d2.push_back(new Node<2>(1, false, 4.0, -3.0));
		nodes2d2.push_back(new Node<2>(2, false, 2.0, -1.0));
		Element<2,2> element2d2(nodes2d2);
		const MatrixDouble *J2d2 = element2d2.GetJacobian();
		TS_ASSERT_DELTA((*J2d2)(0,0), 3.0, 1e-12);
		TS_ASSERT_DELTA((*J2d2)(0,1), 1.0, 1e-12);
		TS_ASSERT_DELTA((*J2d2)(1,0), -1.0, 1e-12);
		TS_ASSERT_DELTA((*J2d2)(1,1), 1.0, 1e-12);
		
		double Det2d = element2d2.GetJacobianDeterminant();
		TS_ASSERT_DELTA(Det2d, 4.0, 1e-12);
		const MatrixDouble *J2d2inv = element2d2.GetInverseJacobian();
		TS_ASSERT_DELTA((*J2d2inv)(0,0), 0.25, 1e-12);
		TS_ASSERT_DELTA((*J2d2inv)(0,1), -0.25, 1e-12);
		TS_ASSERT_DELTA((*J2d2inv)(1,0), 0.25, 1e-12);
		TS_ASSERT_DELTA((*J2d2inv)(1,1), 0.75, 1e-12);
		
		// 3d easy
		std::vector<const Node<3>*> nodes3d;
		nodes3d.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
		nodes3d.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
		nodes3d.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
		nodes3d.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
		Element<3,3> element3d(nodes3d, true);
		const MatrixDouble *J3d = element3d.GetJacobian();
		TS_ASSERT_DELTA((*J3d)(0,0), 1.0, 1e-12);
		TS_ASSERT_DELTA((*J3d)(0,1), 0.0, 1e-12);
		TS_ASSERT_DELTA((*J3d)(0,2), 0.0, 1e-12);
		TS_ASSERT_DELTA((*J3d)(1,0), 0.0, 1e-12);
		TS_ASSERT_DELTA((*J3d)(1,1), 1.0, 1e-12);
		TS_ASSERT_DELTA((*J3d)(1,2), 0.0, 1e-12);
		TS_ASSERT_DELTA((*J3d)(2,0), 0.0, 1e-12);
		TS_ASSERT_DELTA((*J3d)(2,1), 0.0, 1e-12);
		TS_ASSERT_DELTA((*J3d)(2,2), 1.0, 1e-12);
		
		// 3d general
		std::vector<const Node<3>*> nodes3d2;
		nodes3d2.push_back(new Node<3>(0, false, 1.0, 2.0, 3.0));
		nodes3d2.push_back(new Node<3>(1, false, 2.0, 1.0, 3.0));
		nodes3d2.push_back(new Node<3>(2, false, 5.0, 5.0, 5.0));
		nodes3d2.push_back(new Node<3>(3, false, 0.0, 3.0, 4.0));
		Element<3,3> element3d2(nodes3d2, true);
		const MatrixDouble *J3d2 = element3d2.GetJacobian();
		TS_ASSERT_DELTA((*J3d2)(0,0), 1.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2)(0,1), 4.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2)(0,2), -1.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2)(1,0), -1.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2)(1,1), 3.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2)(1,2), 1.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2)(2,0), 0.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2)(2,1), 2.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2)(2,2), 1.0, 1e-4);
		
		double Det3d2 = element3d2.GetJacobianDeterminant();
		TS_ASSERT_DELTA(Det3d2, 7.0, 1e-4);
		const MatrixDouble *J3d2inv = element3d2.GetInverseJacobian();
		TS_ASSERT_DELTA((*J3d2inv)(0,0), 1.0/7.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2inv)(0,1), -6.0/7.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2inv)(0,2), 1.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2inv)(1,0), 1.0/7.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2inv)(1,1), 1.0/7.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2inv)(1,2), 0.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2inv)(2,0), -2.0/7.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2inv)(2,1), -2.0/7.0, 1e-4);
		TS_ASSERT_DELTA((*J3d2inv)(2,2), 1.0, 1e-4);
		
	}
	
    void testNodeToElementConversion(void)
    {
        Point<1> point1(1.0);
        Point<2> point2(2.0,-1.0);
        
        Node<1> node1(0, point1);
        Node<2> node2(0, point2);
        
        Element<0,1> element1(&node1);
        Element<0,2> element2(&node2);
        
        TS_ASSERT_EQUALS(element1.GetNode(0)->GetPoint()[0], point1[0]);
        
        TS_ASSERT_EQUALS(element2.GetNode(0)->GetPoint()[0], point2[0]);
        TS_ASSERT_EQUALS(element2.GetNode(0)->GetPoint()[1], point2[1]);
    }

};

#endif //_TESTELEMENT_HPP_
