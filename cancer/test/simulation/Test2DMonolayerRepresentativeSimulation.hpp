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
#ifndef TEST2DMONOLAYERREPRESENTATIVESIMULATION_HPP_
#define TEST2DMONOLAYERREPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include "TissueSimulation.hpp"
#include "MeshBasedTissue.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "RandomCellKiller.hpp"

/**
 * This class consists of a single test, in which a 2D model
 * of a growing monolayer of cells is simulated for a fixed 
 * period of time.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed. Results can be seen at
 * https://chaste.ediamond.ox.ac.uk/tests.py/profileHistory
 */
class Test2DMonolayerRepresentativeSimulation : public CxxTest::TestSuite
{
public:

    void Test2DMonolayerRepresentativeSimulationForProfiling() throw (Exception)
    {
    	// Set start time
        SimulationTime::Instance()->SetStartTime(0.0);

        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(3.5);

        // Set parameters
        TissueConfig::Instance()->SetMaxTransitGenerations(UINT_MAX);

        // Create some cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (TissueConfig::Instance()->GetStemCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );
            TissueCell cell(TRANSIT, HEALTHY, new StochasticDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Create a force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("Test2DMonolayerRepresentativeSimulationForProfiling");
        simulator.SetEndTime(50.0);

        // Run simulation
		simulator.Solve();

		// Tidy up
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TEST2DMONOLAYERREPRESENTATIVESIMULATION_HPP_*/
