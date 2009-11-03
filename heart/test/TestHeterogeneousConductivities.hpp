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
#ifndef TESTHETEROGENEOUSCONDUCTIVITIES_HPP_
#define TESTHETEROGENEOUSCONDUCTIVITIES_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <fstream>

#include "BidomainProblem.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "Hdf5DataReader.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "SimpleStimulus.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"
#include "PetscSetupAndFinalize.hpp"

using std::ofstream;

/* test class*/
class TestHeterogeneousConductivities : public CxxTest::TestSuite
{
public:
    void TestSimpleSimulation() throw(Exception)
    {
        /*Simulation parameters*/
        HeartConfig::Instance()->SetSimulationDuration(0.7); //ms (falls over after this)
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-6);
        //HeartConfig::Instance()->SetOdeTimeStep(0.01);

        const double width = 0.1;
        const double height = 0.1;
        const double depth = 0.1;

        const unsigned num_elem_x = 8;
        const unsigned num_elem_y = 8;
        const unsigned num_elem_z = 8;

        /* Read the mesh*/
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(num_elem_x, num_elem_y, num_elem_z);
        mesh.Scale(width/num_elem_x, height/num_elem_y, depth/num_elem_z);

        /*Create a cell factory of the type we defined above. */
        GeneralPlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> cell_factory(num_elem_x, width);

        /* monodomain problem class using (a pointer to) the cell factory */
        BidomainProblem<3> problem( &cell_factory );
        problem.SetMesh(&mesh);

        /*tissue properties*/
        std::vector<ChasteCuboid> input_areas;
        std::vector< c_vector<double,3> > intra_conductivities;
        std::vector< c_vector<double,3> > extra_conductivities;
        ChastePoint<3> corner_a(width/2, 0, 0);
        ChastePoint<3> corner_b(width, height, depth);

        input_areas.push_back(ChasteCuboid(corner_a, corner_b));
        //within the cuboid
        intra_conductivities.push_back( Create_c_vector(0.1, 0.1, 0.1) );
        extra_conductivities.push_back( Create_c_vector(0.0, 0.0, 0.0) );
        //This test should *fail* if you comment out the following line
        //(which blocks conductivity on the RHS of the slab).
        HeartConfig::Instance()->SetConductivityHeterogeneities(input_areas, intra_conductivities, extra_conductivities);

        //elsewhere
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.2, 1.2, 1.2));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1.2, 1.2, 1.2));

        /* set  parameters*/
        // HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        // HeartConfig::Instance()->SetCapacitance(1.0);

         /* Output Directory and prefix (for the hdf5 file), relative to CHASTE_TEST_OUTPUT*/
        HeartConfig::Instance()->SetOutputDirectory("slab_results_het_halfcond");
        HeartConfig::Instance()->SetOutputFilenamePrefix("Slab_small");

        /* Initialise the problem*/
        problem.Initialise();

        /* Solve the PDE monodomain equaion*/
        problem.Solve();

        ReplicatableVector voltage_replicated(problem.GetSolution());
        TS_ASSERT_EQUALS(mesh.GetNumNodes() * 2, voltage_replicated.GetSize());
        for (unsigned i=0;i<mesh.GetNumNodes();i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            if (x<width/2)
            {
                 //Left side is stimulated
                 TS_ASSERT_LESS_THAN(-71.0,voltage_replicated[2 * i]);
            }
            else if (x>width/2)
            {
                //Right side is blocked
                TS_ASSERT_LESS_THAN(voltage_replicated[2 * i],-82.0);
            }
        }
    }
};

#endif /*TESTHETEROGENEOUSCONDUCTIVITIES_HPP_*/
