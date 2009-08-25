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
#ifndef TEST2DVERTEXCRYPTREPRESENTATIVESIMULATION_HPP_
#define TEST2DVERTEXCRYPTREPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cancer headers
#include "TissueSimulationArchiver.hpp"

#include "VertexCryptSimulation2d.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "NagaiHondaForce.hpp"
#include "VertexCryptBoundaryForce.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "SloughingCellKiller.hpp"
#include "AbstractCancerTestSuite.hpp"
#include "CancerEventHandler.hpp"

/**
 * This class consists of a single test, in which a 2D vertex model
 * of a colonic crypt is simulated for a fixed period of time.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class Test2DVertexCryptRepresentativeSimulation : public CxxTest::TestSuite
{
public:

    void Test2DVertexCryptRepresentativeSimulationForProfiling() throw (Exception)
    {
        // Create mesh
        unsigned crypt_width = 18;
        unsigned crypt_height = 25;
        Cylindrical2dVertexMesh mesh(crypt_width, crypt_height, true);

        // Create cells
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<mesh.GetNumElements(); elem_index++)
        {
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                 ( TissueConfig::Instance()->GetTransitCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );

            TissueCell cell(TRANSIT, HEALTHY, new SimpleWntCellCycleModel(2));
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> crypt(mesh, cells);

        // Set up Wnt gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(30.0);
        simulator.SetOutputDirectory("Test2DVertexCryptRepresentativeSimulationForProfiling");

        // Make crypt shorter for sloughing
        TissueConfig::Instance()->SetCryptLength(20.0);
        SloughingCellKiller<2> sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Modified timestep to ensure convergence/stability  \todo Make this the default timestep #1098// Modified parameters to make cells equilibriate \todo Make this the default timestep #1098
        simulator.SetDt(0.002);

        // Run simulation
        simulator.Solve();

        // Tidy up
        WntConcentration<2>::Destroy();
    }

};

#endif /*TEST2DMONOLAYERREPRESENTATIVESIMULATION_HPP_*/
