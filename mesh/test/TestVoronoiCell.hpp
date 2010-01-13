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
        vertex1(0) =  -0.2500;
        vertex1(1) = -0.2500;
        vertex1(2) = 1.2500;
        c_vector<double, 3> vertex2;
        vertex2(0) = 1.2500;
        vertex2(1) = -0.2500;
        vertex2(2) = -0.2500;
        c_vector<double, 3> vertex3;
        vertex3(0) = -0.2500;
        vertex3(1) = 1.2500;
        vertex3(2) = -0.2500;
        c_vector<double, 3> vertex4;
        vertex4(0) = 1.2500;
        vertex4(1) = 1.2500;
        vertex4(2) = 1.2500;
        c_vector<double, 3> vertex5;
        vertex5(0) = 1.0;
        vertex5(1) = 1.0;
        vertex5(2) = 1.0;

        Face<3> face1;
        face1.AddVertex(&vertex2);
        face1.AddVertex(&vertex3);
        face1.AddVertex(&vertex4);
        Face<3> face2;
        face2.AddVertex(&vertex1);
        face2.AddVertex(&vertex4);
        face2.AddVertex(&vertex3);
        Face<3> face3;
        face3.AddVertex(&vertex1);
        face3.AddVertex(&vertex2);
        face3.AddVertex(&vertex4);
        Face<3> face4;
        face4.AddVertex(&vertex1);
        face4.AddVertex(&vertex3);
        face4.AddVertex(&vertex2);

        Face<3> face1b;
        face1b.AddVertex(&vertex2);
        face1b.AddVertex(&vertex3);
        face1b.AddVertex(&vertex5);
        Face<3> face2b;
        face2b.AddVertex(&vertex1);
        face2b.AddVertex(&vertex4);
        face2b.AddVertex(&vertex5);
        Face<3> face3b;
        face3b.AddVertex(&vertex1);
        face3b.AddVertex(&vertex2);
        face3b.AddVertex(&vertex5);

        // Face 1 permuted
        Face<3> face1p;
        face1p.AddVertex(&vertex4);
        face1p.AddVertex(&vertex3);
        face1p.AddVertex(&vertex2);

        // Face 1 rotated
        Face<3> face1r;
        face1r.AddVertex(&vertex4);
        face1r.AddVertex(&vertex2);
        face1r.AddVertex(&vertex3);

        VoronoiCell cell1;
        cell1.AddFace(&face1);
        cell1.AddOrientation(true);
        cell1.AddFace(&face2);
        cell1.AddOrientation(true);
        cell1.AddFace(&face3);
        cell1.AddOrientation(true);
        cell1.AddFace(&face4);
        cell1.AddOrientation(true);
        TS_ASSERT_EQUALS(cell1, cell1);

        // A different cell
        VoronoiCell cell1b;
        cell1b.AddFace(&face1b);
        cell1b.AddOrientation(true);
        cell1b.AddFace(&face2b);
        cell1b.AddOrientation(true);
        cell1b.AddFace(&face3b);
        cell1b.AddOrientation(true);
        cell1b.AddFace(&face4);
        cell1b.AddOrientation(true);
        TS_ASSERT_DIFFERS(cell1, cell1b);

        // Like first cell but face 1 permuted
        VoronoiCell cell1p;
        cell1p.AddFace(&face1p);
        cell1p.AddOrientation(true);
        cell1p.AddFace(&face2);
        cell1p.AddOrientation(true);
        cell1p.AddFace(&face3);
        cell1p.AddOrientation(true);
        cell1p.AddFace(&face4);
        cell1p.AddOrientation(true);

        TS_ASSERT_DIFFERS(cell1, cell1p);

        // Like first cell but face 1 rotated, and faces in different order
        VoronoiCell cell1r;
        cell1r.AddFace(&face3);
        cell1r.AddOrientation(true);
        cell1r.AddFace(&face1r);
        cell1r.AddOrientation(true);
        cell1r.AddFace(&face2);
        cell1r.AddOrientation(true);
        cell1r.AddFace(&face4);
        cell1r.AddOrientation(true);

        TS_ASSERT_EQUALS(cell1, cell1r);

        // Null cell
        VoronoiCell cell0;
        TS_ASSERT_DIFFERS(cell1, cell0);
        TS_ASSERT_DIFFERS(cell0, cell1);
        TS_ASSERT_EQUALS(cell0, cell0);

        // Like first cell but face 1 premuted and opposite orientation
        VoronoiCell cell1o;
        cell1o.AddFace(&face3);
        cell1o.AddOrientation(true);
        cell1o.AddFace(&face1p);
        cell1o.AddOrientation(false);
        cell1o.AddFace(&face2);
        cell1o.AddOrientation(true);
        cell1o.AddFace(&face4);
        cell1o.AddOrientation(true);
        TS_ASSERT_EQUALS(cell1, cell1o);
    }

};

#endif /*TESTVORONOICELL_HPP_*/
