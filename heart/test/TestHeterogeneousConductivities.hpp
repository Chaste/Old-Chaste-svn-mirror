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
#ifndef TESTHETEROGENEOUSCONDUCTIVITIES_HPP_
#define TESTHETEROGENEOUSCONDUCTIVITIES_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TrianglesMeshReader.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "Hdf5DataReader.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "SimpleStimulus.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"
#include <iostream>
#include <fstream>
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
       // HeartConfig::Instance()->SetOdeTimeStep(0.01);
   
   
        double num_elem_x=8;
        double num_elem_y=8;
        double num_elem_z=8;
        double width=0.1;
        double height=0.1;
        double depth=0.1;
        
        /* Read the mesh*/
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(num_elem_x, num_elem_y, num_elem_z,  false);
        mesh.Scale(width/num_elem_x, height/num_elem_y, depth/num_elem_z);
   
        
        /*Create a cell factory of the type we defined above. */
        GeneralPlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> cell_factory(num_elem_x, width);
                    
        /* monodomain problem class using (a pointer to) the cell factory */
        MonodomainProblem<3> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);
        
        /*tissue properties*/  
        std::vector< c_vector<double,3> > cornerA;
        std::vector< c_vector<double,3> > cornerB;
        std::vector< c_vector<double,3> > intraConductivities;
        std::vector< c_vector<double,3> > extraConductivities;
        cornerA.push_back( Create_c_vector(width/2, 0, 0) );
        cornerB.push_back( Create_c_vector(width, height, depth) );
        //within the cuboid
        intraConductivities.push_back( Create_c_vector(0.1, 0.1, 0.1) );
        extraConductivities.push_back( Create_c_vector(7.0, 7.0, 7.0) );     
        HeartConfig::Instance()->SetConductivityHeterogeneities(cornerA, cornerB, intraConductivities, extraConductivities); 
        //elsewhere
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.2, 1.2, 1.2));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1.2, 1.2, 1.2));
        
       /* set monodomain parameters*/
       // HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
       // HeartConfig::Instance()->SetCapacitance(1.0);
 
         /* Output Directory and prefix (for the hdf5 file), relative to CHASTE_TEST_OUTPUT*/
        HeartConfig::Instance()->SetOutputDirectory("slab_results_het_halfcond");
        HeartConfig::Instance()->SetOutputFilenamePrefix("Slab_small");    
     
        /*output for MEshalyzer*/
        monodomain_problem.ConvertOutputToMeshalyzerFormat();
        
        /* Initialise the problem*/
        monodomain_problem.Initialise();

        /* Solve the PDE monodomain equaion*/
        monodomain_problem.Solve();
        
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        for (unsigned i=0;i<mesh.GetNumNodes();i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            if (x<=width/2)
            {
                TS_ASSERT_LESS_THAN(-71.0,voltage_replicated[i]);
            }
            else
            {
                TS_ASSERT_LESS_THAN(voltage_replicated[i],-82.0);
            }
        }
        
        
        //Tests of lower-level functionality
        ChastePoint<3> corner1(width/2, 0, 0);
        ChastePoint<3> corner2(width, height, depth);
        ChasteCuboid blocked_region(corner1, corner2); 
        
        ChastePoint<3> p1(0, 0, 0);
        TS_ASSERT_EQUALS(blocked_region.DoesContain(p1), false);
        ChastePoint<3> p2(width/2., height/2., depth/2.);
        TS_ASSERT_EQUALS(blocked_region.DoesContain(p2), true);
        ChastePoint<3> p3(width, height, depth);
        TS_ASSERT_EQUALS(blocked_region.DoesContain(p3), true);
        ChastePoint<2> point2d(width/2., height/2.);
        ChastePoint<1> point1d(width/2.);
        TS_ASSERT_THROWS_ANYTHING(blocked_region.DoesContain(point2d));
        TS_ASSERT_THROWS_ANYTHING(blocked_region.DoesContain(point1d));
        
        
    }
};

#endif /*TESTHETEROGENEOUSCONDUCTIVITIES_HPP_*/
