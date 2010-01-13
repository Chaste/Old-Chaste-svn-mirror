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
#ifndef TESTVORONOIAREAONPERIODICMESH_HPP_
#define TESTVORONOIAREAONPERIODICMESH_HPP_

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>

#include "VoronoiTessellation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TrianglesMeshWriter.hpp"


class TestVoronoiAreaOnPeriodicMesh : public CxxTest::TestSuite
{
public:
    void TestTessellation2NodesOn2dPeriodic() throw (Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 0;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        TS_ASSERT(p_mesh->CheckIsVoronoi());
        TS_ASSERT_DELTA(p_params->GetCryptWidth(),6.0,1e-6);

        // Create Voronoi Tesselation
        VoronoiTessellation<2> tessellation(*p_mesh);

        //  Get two neighbouring nodes on boundary 48 and 53.
        //  Check that they have a common edge
        //  check it is a reasonable length (O(1)?)
        const Face<2> cell_48 = tessellation.rGetFace(48);
        for (unsigned i=0; i<cell_48.GetNumVertices(); i++)
        {
            std::vector< c_vector<double, 2>*> vertices_of_face_48 = cell_48.GetVertices();
            c_vector<double, 2> vertex_of_face_48 = *(vertices_of_face_48[i]);
        }
        const Face<2> cell_53 = tessellation.rGetFace(53);
        for (unsigned i=0; i<cell_53.GetNumVertices(); i++)
        {
            std::vector< c_vector<double, 2>*> vertices_of_face_53 = cell_53.GetVertices();
            c_vector<double, 2> vertex_of_face_53 = *(vertices_of_face_53[i]);
        }

        c_vector<double, 2> location_48 = p_mesh->GetNode(48)->rGetLocation();
        double common_edge_between_48_and_53 = tessellation.GetEdgeLength(48, 53);

        TS_ASSERT_DELTA(tessellation.GetEdgeLength(48, 49), pow(3.0, -0.5), 1e-4);

        TS_ASSERT_DELTA(common_edge_between_48_and_53,  pow(3.0, -0.5), 1e-4);

        //  Check that both cells have a reasonable sized area
        TS_ASSERT_DELTA(tessellation.GetFaceArea(44),  0.5 * pow(3.0, 0.5), 1e-4);
        TS_ASSERT_DELTA(tessellation.GetFacePerimeter(44), 2 * pow(3.0, 0.5), 1e-4);

        TS_ASSERT_DELTA(tessellation.GetFaceArea(48),  0.5 * pow(3.0, 0.5), 1e-4);
        TS_ASSERT_DELTA(tessellation.GetFacePerimeter(48), 2 * pow(3.0, 0.5), 1e-4);
    }
};


#endif /*TESTVORONOIAREAONPERIODICMESH_HPP_*/
