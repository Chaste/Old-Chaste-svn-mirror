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
#ifndef TESTDISTANCEMAPCALCULATOR_
#define TESTDISTANCEMAPCALCULATOR_

#include "TrianglesMeshReader.hpp"
#include "DistanceMapCalculator.hpp"
#include "TetrahedralMesh.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestDistanceMapCalculator : public CxxTest::TestSuite
{
public:
    void TestDistances1D() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");

        std::vector<unsigned> map_origin;
        map_origin.push_back(0u);
        std::vector<double> distances_serial;
        {
            //This is in a block so that we can minimise to scope of the serial mesh (to avoid using it in error)
            TetrahedralMesh<1,1> serial_mesh;
            serial_mesh.ConstructFromMeshReader(mesh_reader);
            TS_ASSERT_EQUALS(serial_mesh.GetNumNodes(), 11u);
            TS_ASSERT_EQUALS(serial_mesh.GetNumElements(), 10u);
            TS_ASSERT_EQUALS(serial_mesh.GetNumBoundaryElements(), 2u);
            DistanceMapCalculator<1,1> distance_calculator_serial(serial_mesh);
            distance_calculator_serial.ComputeDistanceMap(map_origin, distances_serial);
        }
        mesh_reader.Reset();
        
        ParallelTetrahedralMesh<1,1> parallel_mesh;
        parallel_mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS(parallel_mesh.GetNumNodes(), 11u);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumBoundaryElements(), 2u);

        DistanceMapCalculator<1,1> distance_calculator_parallel(parallel_mesh);
        std::vector<double> distances_parallel;
        distance_calculator_parallel.ComputeDistanceMap(map_origin, distances_parallel);
        
        TS_ASSERT_EQUALS(distances_serial.size(), distances_parallel.size());
        for (unsigned index=0; index<distances_parallel.size(); index++)
        {
            try
            {
                c_vector<double, 1> node = parallel_mesh.GetNode(index)->rGetLocation(); //throws if not owned 
                TS_ASSERT_DELTA(distances_serial[index], node(0), 1e-12);
                //TS_ASSERT_DELTA(distances_parallel[index], node(0), 1e-12);
            }
            catch (Exception &e)
            {
            }
        }
    }
        
        
    void TestDistancesToCorner() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/cube_21_nodes_side/Cube21"); // 5x5x5mm cube (internode distance = 0.25mm)

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9261u); // 21x21x21 nodes
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 48000u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4800u);
        
        unsigned far_index=9260u;
        c_vector<double,3> far_corner=mesh.GetNode(far_index)->rGetLocation();
        TS_ASSERT_DELTA( far_corner[0], 0.25, 1e-11);
        TS_ASSERT_DELTA( far_corner[1], 0.25, 1e-11);
        TS_ASSERT_DELTA( far_corner[2], 0.25, 1e-11);
        
        std::vector<unsigned> map_far_corner;
        map_far_corner.push_back(far_index);

        DistanceMapCalculator<3,3> distance_calculator(mesh);
        std::vector<double> distances;
        distance_calculator.ComputeDistanceMap(map_far_corner, distances);

 
        for (unsigned index=0; index<distances.size(); index++)
        {
            c_vector<double, 3> node = mesh.GetNode(index)->rGetLocation();

            double dist = norm_2(far_corner - node);

            TS_ASSERT_DELTA(distances[index], dist, 1e-11);
        }
    }

    void TestDistancesToFace()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/cube_21_nodes_side/Cube21"); // 5x5x5mm cube (internode distance = 0.25mm)

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9261u); // 21x21x21 nodes
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 48000u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4800u);

        DistanceMapCalculator<3,3> distance_calculator(mesh);

        std::vector<unsigned> map_origin;
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            // Get the nodes at the left face of the cube
            if (mesh.GetNode(index)->rGetLocation()[0] + 0.25 < 1e-6)
            {
                map_origin.push_back(index);
            }
        }

        assert(map_origin.size() == 21*21);

        std::vector<double> distances;
        distance_calculator.ComputeDistanceMap(map_origin, distances);

        for (unsigned index=0; index<distances.size(); index++)
        {
            // The distance should be equal to the x-coordinate of the point (minus the offset of the left face of the cube)
            c_vector<double, 3> node = mesh.GetNode(index)->rGetLocation();
            TS_ASSERT_DELTA(distances[index], node[0]+0.25,1e-11);
        }
    }
};

#endif /*TESTDISTANCEMAPCALCULATOR_*/
