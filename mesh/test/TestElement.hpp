#ifndef _TESTELEMENT_HPP_
#define _TESTELEMENT_HPP_

#include "Point.hpp"
#include "Node.hpp"
#include "Element.cpp"
#include "BoundaryElement.cpp"
#include <cxxtest/TestSuite.h>
//#include <iostream>

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
    void TestConstructionForQuadraticBasisFunctions()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.5, 0.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 0.5, 0.0));
        nodes.push_back(new Node<3>(7, false, 0.0, 0.0, 0.5));
        nodes.push_back(new Node<3>(8, false, 0.5, 0.0, 0.5));
        nodes.push_back(new Node<3>(9, false, 0.0, 0.5, 0.5));
        Element<3,3> element(INDEX_IS_NOT_USED, nodes, 2);
        
        // Check nodes on the new element have the right indices
        for (int i=0; i<10; i++)
        {
            TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(i), i);
        }
        
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
    
    void TestConstructionForLinearBasisFunctions()
    {
        std::vector<Node<3>*> corner_nodes;
        corner_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        corner_nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element(INDEX_IS_NOT_USED, corner_nodes, 1);
        
        // Check nodes on the new element have the right indices
        for (int i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(i), i);
        }
        
        for (unsigned i=0; i<corner_nodes.size(); i++)
        {
            delete corner_nodes[i];
        }
    }
    
    
    void TestEquals()
    {
        std::vector<Node<3>*> corner_nodes;
        corner_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        corner_nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element(INDEX_IS_NOT_USED, corner_nodes, 1);
        
        std::vector<Node<3>*> more_nodes;
        more_nodes.push_back(new Node<3>(0, false, 10.0, 10.0, 10.0));
        more_nodes.push_back(new Node<3>(1, false, 11.0, 10.0, 10.0));
        more_nodes.push_back(new Node<3>(2, false, 10.0, 11.0, 10.0));
        more_nodes.push_back(new Node<3>(3, false, 10.0, 10.0, 11.0));
        Element<3,3> another_element(INDEX_IS_NOT_USED, more_nodes, 1);
        
        // test (and cover) equals operator
        another_element = element;
        
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(another_element.GetNode(i)->GetPoint()[j], element.GetNode(i)->GetPoint()[j], 1e-10);
            }
        }
        
        for (unsigned i=0; i<corner_nodes.size(); i++)
        {
            delete corner_nodes[i];
            delete more_nodes[i];
        }
    }
    
    
    void TestGetSetAbstractElementMethods()
    {
        std::vector<Node<3>*> corner_nodes;
        corner_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        corner_nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element(INDEX_IS_NOT_USED, corner_nodes, 1);
        
        TS_ASSERT_EQUALS(element.GetOwnershipSet(),false);
        
        element.SetOwnership(true);
        
        TS_ASSERT_EQUALS(element.GetOwnership(),true);
        TS_ASSERT_EQUALS(element.GetOwnershipSet(),true);
        
        
        for (unsigned i=0; i<corner_nodes.size(); i++)
        {
            delete corner_nodes[i];
        }
        
    }
    
    void TestJacobian()
    {
        // 1d
        std::vector<Node<1>*> nodes1d;
        nodes1d.push_back(new Node<1>(0, false, 2.0));
        nodes1d.push_back(new Node<1>(1, false, 2.5));
        Element<1,1> element1d(INDEX_IS_NOT_USED, nodes1d);
        const c_matrix<double, 1, 1> *J1d = element1d.GetJacobian();
        TS_ASSERT_DELTA((*J1d)(0,0), 0.5, 1e-12);
        
        double Det1d = element1d.GetJacobianDeterminant();
        TS_ASSERT_DELTA(Det1d, 0.5, 1e-12);
        const c_matrix<double, 1, 1> *J1dinv = element1d.GetInverseJacobian();
        TS_ASSERT_DELTA((*J1dinv)(0,0), 2.0, 1e-12);
        
        delete nodes1d[0];
        delete nodes1d[1];
        // 2d easy
        
        std::vector<Node<2>*> nodes2d;
        nodes2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes2d.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes2d.push_back(new Node<2>(2, false, 0.0, 1.0));
        Element<2,2> element2d(INDEX_IS_NOT_USED, nodes2d);
        const c_matrix<double, 2, 2> *J2d = element2d.GetJacobian();
        TS_ASSERT_DELTA((*J2d)(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA((*J2d)(0,1), 0.0, 1e-12);
        TS_ASSERT_DELTA((*J2d)(1,0), 0.0, 1e-12);
        TS_ASSERT_DELTA((*J2d)(1,1), 1.0, 1e-12);
        
        delete nodes2d[0];
        delete nodes2d[1];
        delete nodes2d[2];
        
        //2d general
        std::vector<Node<2>*> nodes2d2;
        nodes2d2.push_back(new Node<2>(0, false, 1.0, -2.0));
        nodes2d2.push_back(new Node<2>(1, false, 4.0, -3.0));
        nodes2d2.push_back(new Node<2>(2, false, 2.0, -1.0));
        Element<2,2> element2d2(INDEX_IS_NOT_USED, nodes2d2);
        const c_matrix<double, 2, 2> *J2d2 = element2d2.GetJacobian();
        TS_ASSERT_DELTA((*J2d2)(0,0), 3.0, 1e-12);
        TS_ASSERT_DELTA((*J2d2)(0,1), 1.0, 1e-12);
        TS_ASSERT_DELTA((*J2d2)(1,0), -1.0, 1e-12);
        TS_ASSERT_DELTA((*J2d2)(1,1), 1.0, 1e-12);
        
        double Det2d = element2d2.GetJacobianDeterminant();
        TS_ASSERT_DELTA(Det2d, 4.0, 1e-12);
        const c_matrix<double, 2, 2> *J2d2inv = element2d2.GetInverseJacobian();
        TS_ASSERT_DELTA((*J2d2inv)(0,0), 0.25, 1e-12);
        TS_ASSERT_DELTA((*J2d2inv)(0,1), -0.25, 1e-12);
        TS_ASSERT_DELTA((*J2d2inv)(1,0), 0.25, 1e-12);
        TS_ASSERT_DELTA((*J2d2inv)(1,1), 0.75, 1e-12);
        
        delete nodes2d2[0];
        delete nodes2d2[1];
        delete nodes2d2[2];
        
        
        // 3d easy
        std::vector<Node<3>*> nodes3d;
        nodes3d.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes3d.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element3d(INDEX_IS_NOT_USED, nodes3d);
        const c_matrix<double, 3, 3> *J3d = element3d.GetJacobian();
        TS_ASSERT_DELTA((*J3d)(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA((*J3d)(0,1), 0.0, 1e-12);
        TS_ASSERT_DELTA((*J3d)(0,2), 0.0, 1e-12);
        TS_ASSERT_DELTA((*J3d)(1,0), 0.0, 1e-12);
        TS_ASSERT_DELTA((*J3d)(1,1), 1.0, 1e-12);
        TS_ASSERT_DELTA((*J3d)(1,2), 0.0, 1e-12);
        TS_ASSERT_DELTA((*J3d)(2,0), 0.0, 1e-12);
        TS_ASSERT_DELTA((*J3d)(2,1), 0.0, 1e-12);
        TS_ASSERT_DELTA((*J3d)(2,2), 1.0, 1e-12);
        
        delete nodes3d[0];
        delete nodes3d[1];
        delete nodes3d[2];
        delete nodes3d[3];
        
        
        // 3d general
        std::vector<Node<3>*> nodes3d2;
        nodes3d2.push_back(new Node<3>(0, false, 1.0, 2.0, 3.0));
        nodes3d2.push_back(new Node<3>(1, false, 2.0, 1.0, 3.0));
        nodes3d2.push_back(new Node<3>(2, false, 5.0, 5.0, 5.0));
        nodes3d2.push_back(new Node<3>(3, false, 0.0, 3.0, 4.0));
        Element<3,3> element3d2(INDEX_IS_NOT_USED, nodes3d2);
        const c_matrix<double, 3, 3> *J3d2 = element3d2.GetJacobian();
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
        const c_matrix<double, 3, 3> *J3d2inv = element3d2.GetInverseJacobian();
        TS_ASSERT_DELTA((*J3d2inv)(0,0), 1.0/7.0, 1e-4);
        TS_ASSERT_DELTA((*J3d2inv)(0,1), -6.0/7.0, 1e-4);
        TS_ASSERT_DELTA((*J3d2inv)(0,2), 1.0, 1e-4);
        TS_ASSERT_DELTA((*J3d2inv)(1,0), 1.0/7.0, 1e-4);
        TS_ASSERT_DELTA((*J3d2inv)(1,1), 1.0/7.0, 1e-4);
        TS_ASSERT_DELTA((*J3d2inv)(1,2), 0.0, 1e-4);
        TS_ASSERT_DELTA((*J3d2inv)(2,0), -2.0/7.0, 1e-4);
        TS_ASSERT_DELTA((*J3d2inv)(2,1), -2.0/7.0, 1e-4);
        TS_ASSERT_DELTA((*J3d2inv)(2,2), 1.0, 1e-4);
        
        delete nodes3d2[0];
        delete nodes3d2[1];
        delete nodes3d2[2];
        delete nodes3d2[3];
    }
    
    void TestNodeToElementConversion(void)
    {
        Point<1> point1(1.0);
        Point<2> point2(2.0,-1.0);
        
        Node<1> node1(0, point1);
        Node<2> node2(0, point2);
        
        BoundaryElement<0,1> element1(INDEX_IS_NOT_USED, &node1);
        BoundaryElement<0,2> element2(INDEX_IS_NOT_USED, &node2);
        
        TS_ASSERT_EQUALS(element1.GetNode(0)->GetPoint()[0], point1[0]);
        
        TS_ASSERT_EQUALS(element2.GetNode(0)->GetPoint()[0], point2[0]);
        TS_ASSERT_EQUALS(element2.GetNode(0)->GetPoint()[1], point2[1]);
    }
    
    void TestGetNodeLocation() throw(Exception)
    {
        Point<2> point1(0.0,1.0);
        Point<2> point2(4.0,6.0);
        Point<2> point3(2.0,3.0);
        
        std::vector<Node<2>*> element_nodes;
        element_nodes.push_back(new Node<2>(0, point1));
        element_nodes.push_back(new Node<2>(0, point2));
        element_nodes.push_back(new Node<2>(0, point3));
        
        Element<2,2> element(INDEX_IS_NOT_USED, element_nodes);
        
        // note that nodes 2 and 3 are swapped by the element constructor
        // to ensure that the jacobian determinant is positive
        TS_ASSERT_EQUALS(element.GetNodeLocation(0,0),0.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(0)(0),0.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(1,0),2.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(1)(0),2.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(2,0),4.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(2)(0),4.0);
        
        delete element_nodes[0];
        delete element_nodes[1];
        delete element_nodes[2];
    }
    
    void TestElementSwapsNodesIfJacobianIsNegative()
    {
        Point<1> a0(0),    a1(1);
        Point<2> b0(0,0),  b1(1,0)  ,b2(0,1);
        Point<3> c0(0,0,0),c1(1,0,0),c2(0,1,0),c3(0,0,1);
        
        Node<1>  na0(0,a0), na1(1,a1);
        Node<2>  nb0(0,b0), nb1(1,b1), nb2(2,b2);
        Node<3>  nc0(0,c0), nc1(1,c1), nc2(2,c2), nc3(3,c3);
        
        
        ////////////////////////////////////////////
        // 1d
        ////////////////////////////////////////////
        std::vector<Node<1>*> nodes_1d_correct;
        nodes_1d_correct.push_back(&na0);
        nodes_1d_correct.push_back(&na1);
        
        std::vector<Node<1>*> nodes_1d_incorrect;
        nodes_1d_incorrect.push_back(&na1);
        nodes_1d_incorrect.push_back(&na0);
        
        Element<1,1>   e_1d_correct_orientation(INDEX_IS_NOT_USED, nodes_1d_correct);
        Element<1,1> e_1d_incorrect_orientation(INDEX_IS_NOT_USED, nodes_1d_incorrect);
        
        // index of second node should be 1
        TS_ASSERT_EQUALS( e_1d_correct_orientation.GetNode(1)->GetIndex(), 1);
        // index of second node for incorrect orientation element should also be 1
        // because the element should have swapped the nodes around
        TS_ASSERT_EQUALS( e_1d_incorrect_orientation.GetNode(1)->GetIndex(), 1);
        
        
        ////////////////////////////////////////////
        // 2d
        ////////////////////////////////////////////
        std::vector<Node<2>*> nodes_2d_correct;
        nodes_2d_correct.push_back(&nb0);
        nodes_2d_correct.push_back(&nb1);
        nodes_2d_correct.push_back(&nb2);
        
        std::vector<Node<2>*> nodes_2d_incorrect;
        nodes_2d_incorrect.push_back(&nb1);
        nodes_2d_incorrect.push_back(&nb0);
        nodes_2d_incorrect.push_back(&nb2);
        
        Element<2,2>   e_2d_correct_orientation(INDEX_IS_NOT_USED, nodes_2d_correct);
        Element<2,2> e_2d_incorrect_orientation(INDEX_IS_NOT_USED, nodes_2d_incorrect);
        
        // index of last node should be 2
        TS_ASSERT_EQUALS( e_2d_correct_orientation.GetNode(2)->GetIndex(), 2);
        
        // index of last node for incorrect orientation element should be 0
        // because the element should have swapped the last two nodes around
        TS_ASSERT_EQUALS( e_2d_incorrect_orientation.GetNode(2)->GetIndex(), 0);
        
        
        ////////////////////////////////////////////
        // 3d
        ////////////////////////////////////////////
        std::vector<Node<3>*> nodes_3d_correct;
        nodes_3d_correct.push_back(&nc0);
        nodes_3d_correct.push_back(&nc1);
        nodes_3d_correct.push_back(&nc2);
        nodes_3d_correct.push_back(&nc3);
        
        std::vector<Node<3>*> nodes_3d_incorrect;
        nodes_3d_incorrect.push_back(&nc0);
        nodes_3d_incorrect.push_back(&nc1);
        nodes_3d_incorrect.push_back(&nc3);
        nodes_3d_incorrect.push_back(&nc2);
        
        Element<3,3>   e_3d_correct_orientation(INDEX_IS_NOT_USED, nodes_3d_correct);
        Element<3,3> e_3d_incorrect_orientation(INDEX_IS_NOT_USED, nodes_3d_incorrect);
        
        // index of last node should be 3
        TS_ASSERT_EQUALS( e_3d_correct_orientation.GetNode(3)->GetIndex(), 3);
        
        // index of last node for incorrect orientation element should be 3
        // because the element should have swapped the last two nodes around
        TS_ASSERT_EQUALS( e_3d_incorrect_orientation.GetNode(3)->GetIndex(), 3);
        
    }
    
    void TestElementCopyConstructor(void)
    {
        std::vector<Node<3>*> cornerNodes;
        cornerNodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        cornerNodes.push_back(new Node<3>(1, false, 1.0, 0.0, 3.0));
        cornerNodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        cornerNodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element(31415, cornerNodes, 1);
        
        // Create a copy of the element and test that it's the same as the original one
        
        Element<3,3> copied_element(element);
        
        TS_ASSERT_EQUALS(element.GetNumNodes(), copied_element.GetNumNodes());
        //check that this version of the copy constructor gives the copied
        //element the same index number
        TS_ASSERT_EQUALS(copied_element.GetIndex(), 31415u);
        
        for (int node = 0; node < 4; node++)
        {
            for (int dimension = 0; dimension < 3; dimension++)
            {
                TS_ASSERT_EQUALS(element.GetNodeLocation(node, dimension), copied_element.GetNodeLocation(node, dimension));
            }
        }
        
        Element<3,3> another_copied_element(element, 2345);
        TS_ASSERT_EQUALS(another_copied_element.GetIndex(), 2345u);
        
        // update a node of the element
        another_copied_element.UpdateNode(1, new Node<3>(4, false, 0.0, 0.0, 2.0));
        TS_ASSERT_EQUALS(another_copied_element.GetNodeLocation(1, 2), 2.0);
    }
    
    void TestBoundaryElement()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.5, 0.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 0.5, 0.0));
        nodes.push_back(new Node<3>(7, false, 0.0, 0.0, 0.5));
        nodes.push_back(new Node<3>(8, false, 0.5, 0.0, 0.5));
        nodes.push_back(new Node<3>(9, false, 0.0, 0.5, 0.5));
        BoundaryElement<3,3> element(INDEX_IS_NOT_USED, nodes, 2);
        
    }
    
    
    void TestCircum1d(void)
    {
        std::vector<Node<1>*> cornerNodes;
        cornerNodes.push_back(new Node<1>(0, false, 10.0));
        cornerNodes.push_back(new Node<1>(1, false, 15.0));
        
        Element<1,1> element(0, cornerNodes, 1);
        
        c_vector <double, 2> circum=element.CalculateCircumsphere();
        TS_ASSERT_DELTA(circum[0], 12.5, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum[1]), 2.5, 1e-7);
        
        TS_ASSERT_DELTA(element.CalculateCircumsphereVolume(), 5.0, 1e-7);
        TS_ASSERT_DELTA(element.CalculateQuality(), 1.0, 1e-7);
    }
    void TestCircum2d(void)
    {
        std::vector<Node<2>*> equilateral_nodes;
        equilateral_nodes.push_back(new Node<2>(0, false, 2.0, 0.0));
        equilateral_nodes.push_back(new Node<2>(1, false, -1.0,sqrt(3)));
        equilateral_nodes.push_back(new Node<2>(2, false, -1.0,-sqrt(3)));
        
        Element<2,2> equilateral_element(0, equilateral_nodes, 1);
        
        c_vector <double, 3> circum=equilateral_element.CalculateCircumsphere();
        TS_ASSERT_DELTA(circum[0], 0.0, 1e-7);
        TS_ASSERT_DELTA(circum[1], 0.0, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum[2]), 2.0, 1e-7);
        
        TS_ASSERT_DELTA(equilateral_element.CalculateCircumsphereVolume(), 4.0*M_PI, 1e-7);
        TS_ASSERT_DELTA(equilateral_element.CalculateQuality(), 1.0, 1e-7);
        
        std::vector<Node<2>*> right_angle_nodes;
        right_angle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        right_angle_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        right_angle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        Element<2,2> right_angle_element(0, right_angle_nodes, 1);
        
        c_vector <double, 3> circum2=right_angle_element.CalculateCircumsphere();
        TS_ASSERT_DELTA(circum2[0], 0.5, 1e-7);
        TS_ASSERT_DELTA(circum2[1], 0.5, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum2[2]), 1.0/sqrt(2.0), 1e-7);
        
        TS_ASSERT_DELTA(right_angle_element.CalculateCircumsphereVolume(), M_PI_2, 1e-7);
        TS_ASSERT_DELTA(right_angle_element.CalculateQuality(), 4.0*sqrt(3.0)/9.0, 1e-7);
        
    }
    
    void TestCircum3d(void)
    {
        std::vector<Node<3>*> cornerNodes;
        cornerNodes.push_back(new Node<3>(0, false,  1.0,  1.0,  1.0));
        cornerNodes.push_back(new Node<3>(1, false, -1.0, -1.0,  1.0));
        cornerNodes.push_back(new Node<3>(2, false, -1.0,  1.0, -1.0));
        cornerNodes.push_back(new Node<3>(3, false,  1.0, -1.0, -1.0));
        
        
        Element<3,3> element(0, cornerNodes, 1);
        
        c_vector <double, 4> circum=element.CalculateCircumsphere();
        TS_ASSERT_DELTA(circum[0], 0.0, 1e-7);
        TS_ASSERT_DELTA(circum[1], 0.0, 1e-7);
        TS_ASSERT_DELTA(circum[2], 0.0, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum[3]), sqrt(3.0), 1e-7);
    
        TS_ASSERT_DELTA(element.CalculateCircumsphereVolume(), 4.0*M_PI*sqrt(3), 1e-7);
        TS_ASSERT_DELTA(element.CalculateCircumsphereVolume(), 4.0*M_PI*sqrt(3), 1e-7);
        TS_ASSERT_DELTA(element.CalculateQuality(), 1.0, 1e-7);
        
        std::vector<Node<3>*> right_angle_nodes;
        right_angle_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        right_angle_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        right_angle_nodes.push_back(new Node<3>(3, false, 0.0, 1.0, 0.0));
        right_angle_nodes.push_back(new Node<3>(2, false, 0.0, 0.0, 1.0));
        
        Element<3,3> right_angle_element(0, right_angle_nodes, 1);
        
        c_vector <double, 4> circum2=right_angle_element.CalculateCircumsphere();
        TS_ASSERT_DELTA(circum2[0], 0.5, 1e-7);
        TS_ASSERT_DELTA(circum2[1], 0.5, 1e-7);
        TS_ASSERT_DELTA(circum2[2], 0.5, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum2[3]), sqrt(3.0)/2.0, 1e-7);
        
        TS_ASSERT_DELTA(right_angle_element.CalculateCircumsphereVolume(), sqrt(3)*M_PI_2, 1e-7);
        TS_ASSERT_DELTA(right_angle_element.GetJacobianDeterminant(), 1.0, 1e-7);
        TS_ASSERT_DELTA(right_angle_element.CalculateQuality(), 0.5, 1e-7);
        
    }
    
    
  
};

#endif //_TESTELEMENT_HPP_
