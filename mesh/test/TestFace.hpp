/*

Copyright (C) University of Oxford, 2005-2010

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
        vertex0(0) =  -0.25*(1+DBL_EPSILON);
        vertex0(1) = -0.2500;
        vertex0(2) = 1.2500;
        c_vector<double, 3> vertex1;
        vertex1(0) =  -0.2500;
        vertex1(1) = -0.2500;
        vertex1(2) = 1.2500;

        // Check that there really is a difference
        TS_ASSERT_DIFFERS(norm_2(vertex0-vertex1), 0.0);
        c_vector<double, 3> vertex2;
        vertex2(0) = 1.2500;
        vertex2(1) = -0.2500;
        vertex2(2) = -0.2500;
        c_vector<double, 3> vertex3;
        vertex3(0) =  -0.2500;
        vertex3(1) = 1.2500;
        vertex3(2) = -0.2500;
        c_vector<double, 3> vertex4;
        vertex4(0) =  1.2500;
        vertex4(1) = 1.2500;
        vertex4(2) = 1.2500;
        c_vector<double, 3> vertex5;
        vertex5(0) =  1.0;
        vertex5(1) = 1.0;
        vertex5(2) = 1.0;
        c_vector<double, 3> vertex6;
        vertex6(0) =  2.0;
        vertex6(1) = 2.0;
        vertex6(2) = 2.0;

        Face<3> face0;
        Face<3> face1;
        face1.AddVertex(&vertex2);
        face1.AddVertex(&vertex3);
        face1.AddVertex(&vertex4);

        Face<3> face2;
        face2.AddVertex(&vertex1);
        face2.AddVertex(&vertex5);
        face2.AddVertex(&vertex6);

        Face<3> face3;
        face3.AddVertex(&vertex0);
        face3.AddVertex(&vertex5);
        face3.AddVertex(&vertex6);

        Face<3> face4;
        face4.AddVertex(&vertex3);
        face4.AddVertex(&vertex4);
        face4.AddVertex(&vertex6);

        Face<3> face5;
        face5.AddVertex(&vertex3);
        face5.AddVertex(&vertex4);
        face5.AddVertex(&vertex2);

        Face<3> face6;
        face6.AddVertex(&vertex4);
        face6.AddVertex(&vertex3);
        face6.AddVertex(&vertex2);

        TS_ASSERT_EQUALS(face1,face1);
        TS_ASSERT_DIFFERS(face1,face2);

        TS_ASSERT_EQUALS(face2, face3);
        TS_ASSERT_DIFFERS(face1, face4);

        TS_ASSERT_EQUALS(face1, face5);
        TS_ASSERT_DIFFERS(face1, face6);
        TS_ASSERT_DIFFERS(face5, face6);

        TS_ASSERT_EQUALS(face0, face0);
        TS_ASSERT_DIFFERS(face1, face0);
        TS_ASSERT_DIFFERS(face0, face1);

        TS_ASSERT_EQUALS(face1, -face6);
    }

    void TestFace2D()
    {
        c_vector<double, 2> vertex0;
        vertex0(0) = 0.0;
        vertex0(1) = 0.0;

        c_vector<double, 2> vertex1;
        vertex1(0) = 1.0;
        vertex1(1) = 0.0;

        c_vector<double, 2> vertex2;
        vertex2(0) = 1.0;
        vertex2(1) = 1.0;

        c_vector<double, 2> vertex3;
        vertex3(0) = 0.0;
        vertex3(1) = 1.0;

        Face<2> face;
        face.AddVertex(&vertex0);
        face.AddVertex(&vertex1);
        face.AddVertex(&vertex2);
        face.AddVertex(&vertex3);

        TS_ASSERT_DELTA(face.GetPerimeter(), 4.0, 1e-12);
        TS_ASSERT_DELTA(face.GetArea(), 1.0, 1e-12);

        // Test the reorder method

        Face<2> reordered_face = -face;
        TS_ASSERT_DIFFERS(reordered_face,face);

        reordered_face.OrderVerticesAntiClockwise();

        TS_ASSERT_EQUALS(reordered_face, face);
    }
};

#endif /*TESTFACE_HPP_*/
