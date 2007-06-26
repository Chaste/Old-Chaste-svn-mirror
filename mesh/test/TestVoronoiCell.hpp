#ifndef TESTVORONOICELL_HPP_
#define TESTVORONOICELL_HPP_

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"

#include <cmath>
#include <vector>

class TestVoronoiCell : public CxxTest::TestSuite
{
public:
    void TestCreateCell()
    {
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
        vertex5(0)= 1.0;      
        vertex5(1)=1.0;
        vertex5(2)=1.0;
        
        Face face1;
        face1.mVertices.push_back(&vertex2);
        face1.mVertices.push_back(&vertex3);
        face1.mVertices.push_back(&vertex4);
        Face face2;
        face2.mVertices.push_back(&vertex1);
        face2.mVertices.push_back(&vertex4);
        face2.mVertices.push_back(&vertex3);
        Face face3;
        face3.mVertices.push_back(&vertex1);
        face3.mVertices.push_back(&vertex2);
        face3.mVertices.push_back(&vertex4);
        Face face4;
        face4.mVertices.push_back(&vertex1);
        face4.mVertices.push_back(&vertex3);
        face4.mVertices.push_back(&vertex2);
                
        Face face1b;
        face1b.mVertices.push_back(&vertex2);
        face1b.mVertices.push_back(&vertex3);
        face1b.mVertices.push_back(&vertex5);
        Face face2b;
        face2b.mVertices.push_back(&vertex1);
        face2b.mVertices.push_back(&vertex4);
        face2b.mVertices.push_back(&vertex5);
        Face face3b;
        face3b.mVertices.push_back(&vertex1);
        face3b.mVertices.push_back(&vertex2);
        face3b.mVertices.push_back(&vertex5);
        // face 1 permuted
        Face face1p;
        face1p.mVertices.push_back(&vertex4);
        face1p.mVertices.push_back(&vertex3);
        face1p.mVertices.push_back(&vertex2);
        // face 1 rotated
        Face face1r;
        face1r.mVertices.push_back(&vertex4);
        face1r.mVertices.push_back(&vertex2);
        face1r.mVertices.push_back(&vertex3);
        
        VoronoiCell cell1;
        cell1.mFaces.push_back(&face1);
        cell1.mOrientations.push_back(true);
        cell1.mFaces.push_back(&face2);
        cell1.mOrientations.push_back(true);
        cell1.mFaces.push_back(&face3);
        cell1.mOrientations.push_back(true);
        cell1.mFaces.push_back(&face4);
        cell1.mOrientations.push_back(true);
        TS_ASSERT_EQUALS(cell1, cell1);
        
       // a different cell
        VoronoiCell cell1b;
        cell1b.mFaces.push_back(&face1b);
        cell1b.mOrientations.push_back(true);
        cell1b.mFaces.push_back(&face2b);
        cell1b.mOrientations.push_back(true);
        cell1b.mFaces.push_back(&face3b);
        cell1b.mOrientations.push_back(true);
        cell1b.mFaces.push_back(&face4);
        cell1b.mOrientations.push_back(true);
        TS_ASSERT_DIFFERS(cell1, cell1b);
        
        // like first cell but face 1 permuted
        VoronoiCell cell1p;
        cell1p.mFaces.push_back(&face1p);
        cell1p.mOrientations.push_back(true);
        cell1p.mFaces.push_back(&face2);
        cell1p.mOrientations.push_back(true);
        cell1p.mFaces.push_back(&face3);
        cell1p.mOrientations.push_back(true);
        cell1p.mFaces.push_back(&face4);
        cell1p.mOrientations.push_back(true);
        
        TS_ASSERT_DIFFERS(cell1, cell1p);
        
        // like first cell but face 1 rotated, and faces in different order
        VoronoiCell cell1r;
        cell1r.mFaces.push_back(&face3);
        cell1r.mOrientations.push_back(true);
        cell1r.mFaces.push_back(&face1r);
        cell1r.mOrientations.push_back(true);
        cell1r.mFaces.push_back(&face2);
        cell1r.mOrientations.push_back(true);
        cell1r.mFaces.push_back(&face4);
        cell1r.mOrientations.push_back(true);
        
        TS_ASSERT_EQUALS(cell1, cell1r);
        
        // null cell
        VoronoiCell cell0;
        TS_ASSERT_DIFFERS(cell1, cell0);
        TS_ASSERT_DIFFERS(cell0, cell1);
        TS_ASSERT_EQUALS(cell0, cell0);
        
        // like first cell but face 1 premuted and opposite orientation
        VoronoiCell cell1o;
        cell1o.mFaces.push_back(&face3);
        cell1o.mOrientations.push_back(true);
        cell1o.mFaces.push_back(&face1p);
        cell1o.mOrientations.push_back(false);
        cell1o.mFaces.push_back(&face2);
        cell1o.mOrientations.push_back(true);
        cell1o.mFaces.push_back(&face4);
        cell1o.mOrientations.push_back(true);
        TS_ASSERT_EQUALS(cell1, cell1o);
        
    }
    
};


#endif /*TESTVORONOICELL_HPP_*/
