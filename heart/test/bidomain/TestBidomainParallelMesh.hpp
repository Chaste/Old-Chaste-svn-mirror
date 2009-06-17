/*

Copyright (C) University of Oxford, 2005-2009

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

#ifndef TESTBIDOMAINPARALLELMESH_HPP_
#define TESTBIDOMAINPARALLELMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BidomainProblem.hpp"
#include "DistributedVector.hpp"
#include "HeartConfig.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestBidomainParallelMesh : public CxxTest::TestSuite
{
public:

    void TestBidomainProblemWithDistributedMesh2D()
    {
        HeartConfig::Instance()->SetSimulationDuration(1);  //ms
        HeartConfig::Instance()->SetOutputDirectory("DistributedMesh2d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("tetrahedral2d");

        // The default stimulus in PlaneStimulusCellFactory is not enough to generate propagation
        // here, increasing it an order of magnitude
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-6000);

        // To avoid an issue with the Event handler only one simulation should be
        // in existance at a time: therefore monodomain simulation is defined in a block
        double seq_ave_voltage;
        {
            ///////////////////////////////////////////////////////////////////
            // TetrahedralMesh
            ///////////////////////////////////////////////////////////////////
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            BidomainProblem<2> nondistributed_problem( &cell_factory );
            nondistributed_problem.SetMesh(&mesh);
            nondistributed_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            nondistributed_problem.Solve();

            DistributedVector dist_nondistributed_voltage = nondistributed_problem.GetSolutionDistributedVector();
            DistributedVector::Stripe nondistributed_voltage(dist_nondistributed_voltage, 0);
            DistributedVector::Stripe nondistributed_potential(dist_nondistributed_voltage, 1);

            double seq_local_ave_voltage = 0.0;

            for (DistributedVector::Iterator index = dist_nondistributed_voltage.Begin();
                 index != dist_nondistributed_voltage.End();
                 ++index)
            {
                if (index.Global==0)
                {
                    TS_ASSERT_LESS_THAN(0, nondistributed_voltage[index]);
                }

                seq_local_ave_voltage += nondistributed_voltage[index];
            }

            MPI_Reduce(&seq_local_ave_voltage, &seq_ave_voltage, 1, MPI_DOUBLE, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);
            seq_ave_voltage /= mesh.GetNumNodes();
        }


        ///////////////////////////////////////////////////////////////////
        // ParallelTetrahedralMesh
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputFilenamePrefix("distributed1d");

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        ParallelTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        BidomainProblem<2> distributed_problem( &cell_factory );

        distributed_problem.SetMesh(&mesh);

        distributed_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        distributed_problem.Solve();

        DistributedVector dist_distributed_voltage = distributed_problem.GetSolutionDistributedVector();
        DistributedVector::Stripe distributed_voltage(dist_distributed_voltage, 0);
        DistributedVector::Stripe distributed_potential(dist_distributed_voltage, 1);

        double para_local_ave_voltage = 0.0;

        for (DistributedVector::Iterator index = dist_distributed_voltage.Begin();
             index != dist_distributed_voltage.End();
             ++index)
        {
            if (index.Global==0)
            {
                TS_ASSERT_LESS_THAN(0, distributed_voltage[index]);
            }

            para_local_ave_voltage += distributed_voltage[index];
        }


        double para_ave_voltage;
        MPI_Reduce(&para_local_ave_voltage, &para_ave_voltage, 1, MPI_DOUBLE, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);
        para_ave_voltage /= mesh.GetNumNodes();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        if (PetscTools::AmMaster())
        {
            std::cout << seq_ave_voltage << "  " << para_ave_voltage << std::endl;
            TS_ASSERT_DELTA(seq_ave_voltage, para_ave_voltage, 1.0);
        }

    }
};

#endif /*TESTBIDOMAINPARALLELMESH_HPP_*/
