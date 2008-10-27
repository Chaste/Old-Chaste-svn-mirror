/*

Copyright (C) University of Oxford, 2008

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

#include "VoronoiCell.hpp"
#include "VoronoiTessellation.hpp"
#include "TetrahedralMesh.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "Exception.hpp"
#include "TrianglesMeshWriter.hpp"


class TestVoronoiAreaOnPeriodicMesh : public CxxTest::TestSuite
{
public:
    void TestTessellation2NodesOn2dPeriodic() throw (Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 0;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        TS_ASSERT(p_mesh->CheckVoronoi());
        TS_ASSERT_DELTA(p_params->GetCryptWidth(),6.0,1e-6);

        // Create Voronoi Tesselation
        VoronoiTessellation<2> tessellation(*p_mesh);

        //  Get two neighbouring nodes on boundary 48 and 53.
        //  Check that they have a common edge
        //  check it is a reasonable length (O(1)?)
        const Face<2> cell48 = *(tessellation.GetFace(48u));
        for (unsigned i=0; i<cell48.GetNumVertices(); i++)
        {
            std::vector< c_vector<double, 2>*> vertices_of_face48 = cell48.GetVertices();
            c_vector<double, 2> vertex_of_face48 = *(vertices_of_face48[i]);
        }
        const Face<2> cell53 = *(tessellation.GetFace(53u));
        for (unsigned i=0; i<cell53.GetNumVertices(); i++)
        {
            std::vector< c_vector<double, 2>*> vertices_of_face53 = cell53.GetVertices();
            c_vector<double, 2> vertex_of_face53 = *(vertices_of_face53[i]);
        }

        c_vector<double, 2> location48 = p_mesh->GetNode(48)->rGetLocation();
        double common_edge_between48and53 = tessellation.GetEdgeLength(48u, 53u);

        TS_ASSERT_DELTA(tessellation.GetEdgeLength(48u, 49u), pow(3.0, -0.5), 1e-4);

        TS_ASSERT_DELTA(common_edge_between48and53,  pow(3.0, -0.5), 1e-4);

        //  Check that both cells have a reasonable sized area
        TS_ASSERT_DELTA(tessellation.GetFaceArea(44u),  0.5 * pow(3.0, 0.5), 1e-4);
        TS_ASSERT_DELTA(tessellation.GetFacePerimeter(44u), 2 * pow(3.0, 0.5), 1e-4);

        TS_ASSERT_DELTA(tessellation.GetFaceArea(48u),  0.5 * pow(3.0, 0.5), 1e-4);
        TS_ASSERT_DELTA(tessellation.GetFacePerimeter(48u), 2 * pow(3.0, 0.5), 1e-4);
    }
};


#endif /*TESTVORONOIAREAONPERIODICMESH_HPP_*/
