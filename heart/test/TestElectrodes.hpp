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

#ifndef TESTELECTRODES_HPP_
#define TESTELECTRODES_HPP_


#include <cxxtest/TestSuite.h>
#include <vector>
#include "Electrodes.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedVector.hpp"

class TestElectrodes : public CxxTest::TestSuite
{
public:
    void TestElectrodeGrounded2dAndSwitchOff() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ParallelTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double magnitude = 543.324;
        double duration = 2.0; //ms
        Electrodes<2> electrodes(mesh,true,0,0.0,10.0,magnitude,duration);

        BoundaryConditionsContainer<2,2,2>* p_bcc = electrodes.GetBoundaryConditionsContainer();

        for(ParallelTetrahedralMesh<2,2>::BoundaryElementIterator iter
           = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[0] - 0.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);
                TS_ASSERT_DELTA(value,magnitude,1e-12);
            }
        }

        unsigned num_grounded_nodes = 0u;

        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter=mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            Node<2>* p_node = &(*iter);
            if (p_bcc->HasDirichletBoundaryCondition(p_node, 1))
            {
                double x_val = p_node->rGetLocation()[0];
                TS_ASSERT_DELTA(x_val, 10.0, 1e-12);
                num_grounded_nodes++;
                TS_ASSERT_EQUALS(p_bcc->GetDirichletBCValue(p_node, 1), 0.0);
            }
        }

        unsigned num_grounded_nodes_reduced;
        int mpi_ret = MPI_Allreduce(&num_grounded_nodes, &num_grounded_nodes_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        TS_ASSERT_EQUALS(mpi_ret, MPI_SUCCESS);

        TS_ASSERT_EQUALS(num_grounded_nodes_reduced, 11u);

        TS_ASSERT_THROWS_THIS(Electrodes<2> bad_electrodes(mesh,true,0,5.0,10.0,magnitude,duration),
                "Minimum value of coordinate is not the value given");
        TS_ASSERT_THROWS_THIS(Electrodes<2> bad_electrodes(mesh,true,0,0.0,30.0,magnitude,duration),
                "Maximum value of coordinate is not the value given");

        TS_ASSERT_EQUALS(electrodes.SwitchOff(0.0), false); // t<end time
        TS_ASSERT_EQUALS(electrodes.SwitchOff(1.0), false); // t<end time
        TS_ASSERT_EQUALS(electrodes.SwitchOff(1.9), false); // t<end time
        TS_ASSERT_EQUALS(electrodes.SwitchOff(1.99),false); // t<end time
        TS_ASSERT_EQUALS(electrodes.SwitchOff(2.0+1e-12), true); // true as t>end_time
        TS_ASSERT_EQUALS(electrodes.SwitchOff(2.1), false); // false as electrodes has been switched off
        TS_ASSERT_EQUALS(electrodes.SwitchOff(4.0), false); // false as electrodes has been switched off
    }


    void TestElectrodeUngrounded2d() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ParallelTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double magnitude = 543.324;
        double duration = 2.0;
        Electrodes<2> electrodes(mesh,false,0,0,10,magnitude,duration);

        BoundaryConditionsContainer<2,2,2>* p_bcc = electrodes.GetBoundaryConditionsContainer();

        for(ParallelTetrahedralMesh<2,2>::BoundaryElementIterator iter
                = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[0] - 0.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);

                TS_ASSERT_DELTA(value,magnitude,1e-12);
            }


            if ( fabs((*iter)->CalculateCentroid()[0] - 10.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);
                TS_ASSERT_DELTA(value,-magnitude,1e-12);
            }
        }
    }

    void TestElectrodeGrounded3d() throw (Exception)
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(10,10,10);

        double magnitude = 543.324;
        double duration = 2.0;
        Electrodes<3> electrodes(mesh,true,1,0,10,magnitude,duration);

        BoundaryConditionsContainer<3,3,2>* p_bcc = electrodes.GetBoundaryConditionsContainer();

        for(ParallelTetrahedralMesh<3,3>::BoundaryElementIterator iter
                = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[1] - 0.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);

                TS_ASSERT_DELTA(value,magnitude,1e-12);
            }
        }

        unsigned num_grounded_nodes = 0u;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            Node<3>* p_node = mesh.GetNode(i);
            if (p_bcc->HasDirichletBoundaryCondition(p_node, 1))
            {
                double y_val = p_node->rGetLocation()[1];
                TS_ASSERT_DELTA(y_val, 10.0, 1e-12);
                num_grounded_nodes++;
                TS_ASSERT_EQUALS(p_bcc->GetDirichletBCValue(p_node, 1), 0.0);
            }
        }
        TS_ASSERT_EQUALS(num_grounded_nodes, 121u);
    }

    void TestElectrodeUngrounded3d() throw (Exception)
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(10,10,10);

        double magnitude = 543.324;
        double duration = 2.0;
        Electrodes<3> electrodes(mesh,false,1,0,10,magnitude,duration);

        BoundaryConditionsContainer<3,3,2>* p_bcc = electrodes.GetBoundaryConditionsContainer();

        for(ParallelTetrahedralMesh<3,3>::BoundaryElementIterator iter
                = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[1] - 0.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);

                TS_ASSERT_DELTA(value,magnitude,1e-12);
            }


            if ( fabs((*iter)->CalculateCentroid()[1] - 10.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);
                TS_ASSERT_DELTA(value,-magnitude,1e-12);
            }
        }
    }
};



#endif /*TESTELECTRODES_HPP_*/
