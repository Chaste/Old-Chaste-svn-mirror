#ifndef TESTFACE_HPP_
#define TESTFACE_HPP_

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "Face.hpp"


class TestFace : public CxxTest::TestSuite
{
public:
    void TestEquality()
    {
        c_vector<double, 3> vertex0;
        vertex0(0)= -0.25000000000000001;      
        vertex0(1)=-0.2500;
        vertex0(2)=1.2500;
        c_vector<double, 3> vertex1;
        vertex1(0)= -0.2500;      
        vertex1(1)=-0.2500;
        vertex1(2)=1.2500;
        c_vector<double, 3> vertex2;        
        vertex2(0)=1.2500;
        vertex2(1)=-0.2500;
        vertex2(2)=-0.2500;
        c_vector<double, 3> vertex3; 
        vertex3(0)= -0.2500;      
        vertex3(1)=1.2500;
        vertex3(2)=-0.2500;
        c_vector<double, 3> vertex4;
        vertex4(0)= 1.2500;      
        vertex4(1)=1.2500;
        vertex4(2)=1.2500;
        c_vector<double, 3> vertex5;
        vertex4(0)= 1.0;      
        vertex4(1)=1.0;
        vertex4(2)=1.0;
        c_vector<double, 3> vertex6;
        vertex4(0)= 2.0;      
        vertex4(1)=2.0;
        vertex4(2)=2.0;
        Face face0;
        Face face1;
        face1.mVertices.push_back(&vertex2);
        face1.mVertices.push_back(&vertex3);
        face1.mVertices.push_back(&vertex4);
        Face face2;
        face2.mVertices.push_back(&vertex1);
        face2.mVertices.push_back(&vertex5);
        face2.mVertices.push_back(&vertex6);
        Face face3;
        face3.mVertices.push_back(&vertex0);
        face3.mVertices.push_back(&vertex5);
        face3.mVertices.push_back(&vertex6);
        Face face4;
        face4.mVertices.push_back(&vertex3);
        face4.mVertices.push_back(&vertex4);
        face4.mVertices.push_back(&vertex6);
        Face face5;
        face5.mVertices.push_back(&vertex3);
        face5.mVertices.push_back(&vertex4);
        face5.mVertices.push_back(&vertex2);
        Face face6;
        face6.mVertices.push_back(&vertex4);
        face6.mVertices.push_back(&vertex3);
        face6.mVertices.push_back(&vertex2);
        
        TS_ASSERT_EQUALS(face1,face1);
        TS_ASSERT_DIFFERS(face1,face2);
        
        TS_ASSERT_THROWS_ANYTHING(face2 == face3); //Bad bug
        TS_ASSERT_THROWS_ANYTHING(face1 != face4); //Bad bug
        
        TS_ASSERT_EQUALS(face1, face5);
        TS_ASSERT_DIFFERS(face1, face6);
        TS_ASSERT_DIFFERS(face5, face6);
        
        TS_ASSERT_EQUALS(face0, face0);
        TS_ASSERT_DIFFERS(face1, face0);
        TS_ASSERT_DIFFERS(face0, face1);
        
        TS_ASSERT_EQUALS(face1, -face6);
    }
};

#endif /*TESTFACE_HPP_*/
