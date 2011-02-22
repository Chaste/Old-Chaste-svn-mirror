/*

Copyright (C) University of Oxford, 2005-2011

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





/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTBIDOMAINWITHBATHTUTORIAL_HPP_
#define TESTBIDOMAINWITHBATHTUTORIAL_HPP_
/*
 * = An example showing how to run a bidomain simulation for tissue contained in a perfusing bath =
 * 
 * In this tutorial we show how the changes the need to be made when running a simulation of 
 * cardiac tissue contained in a bath.
 *
 */

/* The usual headers are included */
#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"
/* This test will show how to load a mesh in the test and pass it into the problem,
 * for which the following includes are needed */
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

/* Define the test */
class TestRunningBidomainSimulationsTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!
    
    void TestWithBathAndElectrodes() throw (Exception)
    {
        /* First, set the end time and output info. In this simulation
         * we'll explicitly read the mesh, alter it, then pass it
         * to the problem class, so we don't set the mesh file name.
         */
        HeartConfig::Instance()->SetSimulationDuration(3.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainTutorialWithBath");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        /* Bath problems seem to require decreased ODE timesteps.
         */
        HeartConfig::Instance()->SetOdeTimeStep(0.001);  //ms

        /* Use the {{{PlaneStimulusCellFactory}}} to define a set
         * of Luo-Rudy cells. We pass the stimulus magnitude as 0.0 
         * as we don't want any stimulated cells
         */
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML,2> cell_factory(0.0);

        /*
         * Now, we load up a rectangular mesh (in triangle/tetgen format), done as follows,
         * using {{{TrianglesMeshReader}}}.
         */
        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        /* In bath problems, each element has an attribute which must be set
         * to 0 (cardiac tissue) or 1 (bath). This can be done by having an
         * extra column in the element file (see the file formats documentation, 
         * or for example 
         * mesh/test/data/1D_0_to_1_10_elements_with_two_attributes.ele,
         * and note that the header in this file has 1 at the end to indicate that
         * the file defines an attribute for each element. We have read in a mesh
         * without this type of information set up, so we set it up manually,
         * by looping over elements and setting those more than 2mm from the centre
         * as bath elements (by default, the others are cardiac elements).
         */
        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            double y = mesh.GetElement(i)->CalculateCentroid()[1];
            if( sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05)) > 0.02 )
            {
                mesh.GetElement(i)->SetRegion(HeartRegionCode::BathRegion());
            }
        }

        /* Now we define the electrodes. First define the magnitude of the electrodes
         * (ie the magnitude of the boundary extracellular stimulus), and the duration
         * it lasts for. Currently, electrodes switch on at time 0 and have constant magnitude
         * until they are switched off. (Note that this test has a small range of
         * magnitudes that will work, perhaps because the electrodes are close to the tissue).
         */
        //-1e4 is under threshold, -1.4e4 too high - crashes the cell model
        double magnitude = -1.1e4; // uA/cm^2
        double start_time = 0.0;
        double duration = 2; //ms

        /* Electrodes work in two ways: the first electrode applies an input flux, and
         * the opposite electrode can either be grounded or apply an equal and opposite
         * flux (ie an output flux). The `false` here indicates the second electrode
         * is not grounded, ie has an equal and opposite flux. The "0" indicates
         * that the electrodes should be applied to the bounding surfaces in the x-direction
         * (1 would be y-direction, 2 z-direction), which are X=0.0 and X=0.1 in the given mesh.
         * (This explains why the full mesh ought to be rectangular/cuboid - the nodes on 
         * x=xmin and x=xmax ought to be form two surfaces of equal area.
         */
        HeartConfig::Instance()->SetElectrodeParameters(false, 0, magnitude, start_time, duration);

        /* Now create the problem class, using the cell factory and passing
         * in `true` as the second argument to indicate we are solving a bath
         * problem..
         */
        BidomainProblem<2> bidomain_problem( &cell_factory, true );

        /* ..set the mesh and electrodes.. */
        bidomain_problem.SetMesh(&mesh);

        /* ..and solve as before. */
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        /* The results can be visualised as before. '''Note:''' The voltage is only
         * defined at cardiac nodes (a node contained in ''any'' cardiac element), but
         * for visualisation and computation a 'fake' value of ZERO is given for the
         * voltage at bath nodes.
         *
         * EMPTYLINE
         *
         * Finally, we can check that an AP was induced in any of the cardiac
         * cells. We use a `ReplicatableVector` as before, and make sure we
         * only check the voltage at cardiac cells.
         */
        Vec solution = bidomain_problem.GetSolution(); // the Vs and phi_e's, as a PetSc vector
        ReplicatableVector solution_repl(solution);

        bool ap_triggered = false;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (HeartRegionCode::IsRegionTissue( mesh.GetNode(i)->GetRegion() ))
            {
                if (solution_repl[2*i] > 0.0) // 2*i, ie the voltage for this node (would be 2*i+1 for phi_e for this node)
                {
                    ap_triggered = true;
                }
            }
        }
        TS_ASSERT(ap_triggered);
    }
};

#endif /*TESTBIDOMAINWITHBATHTUTORIAL_HPP_*/
